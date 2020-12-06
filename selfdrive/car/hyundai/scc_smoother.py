import copy
from random import randint
import numpy as np
from common.numpy_fast import clip, interp
from cereal import car, log
from selfdrive.car.hyundai.interface import ButtonType
from selfdrive.config import Conversions as CV, RADAR_TO_CAMERA
from selfdrive.car.hyundai.values import Buttons
from common.params import Params
from selfdrive.controls.lib.drive_helpers import V_CRUISE_MAX, V_CRUISE_MIN, V_CRUISE_DELTA_KM, V_CRUISE_DELTA_MI
from selfdrive.road_speed_limiter import road_speed_limiter_get_max_speed

# do not modify
MIN_SET_SPEED = V_CRUISE_MIN
MAX_SET_SPEED = V_CRUISE_MAX

ALIVE_COUNT_MIN = 6
ALIVE_COUNT_MAX = 8

WAIT_COUNT_MIN = 8
WAIT_COUNT_MAX = 15

MAX_ACC_BUFFER_COUNT = 5
LIMIT_ACCEL = 10.
LIMIT_DECEL = 18.

EventName = car.CarEvent.EventName

ButtonType = car.CarState.ButtonEvent.Type
ButtonPrev = ButtonType.unknown
ButtonCnt = 0
LongPressed = False

class CruiseState:
  STOCK = 0
  SMOOTH = 1
  COUNT = 2

class SccSmoother:
  def __init__(self, accel_gain, decel_gain, curvature_gain):

    self.accel_gain = accel_gain
    self.decel_gain = decel_gain
    self.curvature_gain = curvature_gain

    self.last_cruise_buttons = Buttons.NONE
    self.target_speed = 0

    self.started_frame = 0
    self.accel_buf = []
    self.max_set_speed_buf = []
    self.max_set_speed = 0
    self.wait_timer = 0
    self.alive_timer = 0
    self.btn = Buttons.NONE

    self.alive_count = ALIVE_COUNT_MIN
    self.wait_count = WAIT_COUNT_MIN

    self.path_x = np.arange(192)
    self.curve_speed = 0.

    self.state_changed_for_events = False

    self.state = int(Params().get('SccSmootherState'))
    self.scc_smoother_enabled = Params().get('SccSmootherEnabled') == b'1'
    self.slow_on_curves = Params().get('SccSmootherSlowOnCurves') == b'1'

    self.sync_set_speed_while_gas_pressed = True

  def reset(self):
    self.accel_buf = []
    self.max_set_speed_buf = []
    self.max_set_speed = 0
    self.wait_timer = 0
    self.alive_timer = 0
    self.btn = Buttons.NONE
    self.target_speed = 0

  def create_clu11(self, packer, frame, bus, clu11, button):
    values = copy.copy(clu11)
    values["CF_Clu_CruiseSwState"] = button
    values["CF_Clu_AliveCnt1"] = frame
    return packer.make_can_msg("CLU11", bus, values)

  def is_active(self, frame):
    return frame - self.started_frame <= self.alive_count + self.wait_count

  def dispatch_cancel_buttons(self, CC, CS):
    changed = False
    if self.last_cruise_buttons != CS.cruise_buttons:
      self.last_cruise_buttons = CS.cruise_buttons

      if not CS.cruiseState_enabled:
        if CS.cruise_buttons == Buttons.CANCEL:
          self.state += 1
          if self.state >= CruiseState.COUNT:
            self.state = 0

          Params().put('SccSmootherState', str(self.state))
          self.state_changed_for_events = True
          changed = True

    CC.sccSmoother.state = self.state
    return changed

  def inject_events(self, events):
    if self.state_changed_for_events:
      self.state_changed_for_events = False
      events.add(EventName.sccSmootherStatus)

  def cal_max_speed(self, frame, CC, CS, sm, clu11_speed):

    # Limits on curvature and road speed limit

    if self.max_set_speed == 0 or frame % 5 == 0:

      self.curve_speed = self.get_curve_speed(sm, clu11_speed * CV.KPH_TO_MS) * CV.MS_TO_KPH

      max_speed = min(CC.cruiseOpMaxSpeed, self.curve_speed)
      max_speed = min(max_speed, road_speed_limiter_get_max_speed(CS, max_speed))

      self.max_set_speed_buf.append(max_speed)
      if len(self.max_set_speed_buf) > 5:
        self.max_set_speed_buf.pop(0)

      self.max_set_speed = sum(self.max_set_speed_buf) / len(self.max_set_speed_buf)

  def update(self, enabled, can_sends, packer, CC, CS, frame, apply_accel, controls):

    if not self.scc_smoother_enabled:
      return

    clu11_speed = CS.clu11["CF_Clu_Vanz"]

    if self.dispatch_cancel_buttons(CC, CS):
      return

    if self.state == CruiseState.STOCK or not CS.acc_mode or \
        not enabled or not CS.cruiseState_enabled or CS.cruiseState_speed < 1. or \
        CS.cruiseState_speed > 254 or CS.standstill or \
        CS.cruise_buttons != Buttons.NONE or \
        CS.brake_pressed:

      CC.sccSmoother.logMessage = '{:.2f},{:d},{:d},{:d},{:d},{:.1f},{:d},{:d},{:d}' \
        .format(float(apply_accel*CV.MS_TO_KPH), int(CS.acc_mode), int(enabled), int(CS.cruiseState_enabled), int(CS.standstill), float(CS.cruiseState_speed),
                int(CS.cruise_buttons), int(CS.brake_pressed), int(CS.gas_pressed))

      self.reset()
      self.wait_timer = ALIVE_COUNT_MAX + WAIT_COUNT_MAX
      return

    current_set_speed = CS.cruiseState_speed * CV.MS_TO_KPH

    accel, override_acc = self.cal_acc(apply_accel, CS, clu11_speed, controls.sm)

    if CS.gas_pressed:
      self.target_speed = clu11_speed
      if clu11_speed > controls.cruiseOpMaxSpeed and self.sync_set_speed_while_gas_pressed:
        set_speed = clip(clu11_speed, MIN_SET_SPEED, MAX_SET_SPEED)
        CC.cruiseOpMaxSpeed = controls.cruiseOpMaxSpeed = controls.v_cruise_kph = set_speed
    else:
      self.target_speed = clu11_speed + accel

    self.cal_max_speed(frame, CC, CS, controls.sm, clu11_speed)
    self.target_speed = clip(self.target_speed, MIN_SET_SPEED, self.max_set_speed)

    CC.sccSmoother.logMessage = '{:.2f}/{:.2f}, {:.2f}, {:.1f}/{:d}, btn:{:d}' \
      .format(float(apply_accel * CV.MS_TO_KPH), float(override_acc), float(accel), float(self.target_speed),
              int(self.curve_speed), int(self.btn))

    #CC.sccSmoother.logMessage = '{:.2f}/{:.2f}, {:.2f}, {:.1f}/{:d}, btn:{:d}' \
    #  .format(float(apply_accel*CV.MS_TO_KPH), float(override_acc), float(accel), float(self.target_speed), int(current_set_speed), int(self.btn))

    if self.wait_timer > 0:
      self.wait_timer -= 1
    else:

      if self.alive_timer == 0:
        self.btn = self.get_button(clu11_speed, current_set_speed)
        self.alive_count = randint(ALIVE_COUNT_MIN, ALIVE_COUNT_MAX)

      if self.btn != Buttons.NONE:
        can_sends.append(self.create_clu11(packer, self.alive_timer, CS.scc_bus, CS.clu11, self.btn))

        if self.alive_timer == 0:
          self.started_frame = frame

        self.alive_timer += 1

        if self.alive_timer >= self.alive_count:
          self.alive_timer = 0
          self.wait_count = randint(WAIT_COUNT_MIN, WAIT_COUNT_MAX)
          self.wait_timer = self.wait_count
          self.btn = Buttons.NONE


  def get_button(self, clu11_speed, current_set_speed):

    error = self.target_speed - current_set_speed

    if abs(error) < 0.9:
      return Buttons.NONE

    return Buttons.RES_ACCEL if error > 0 else Buttons.SET_DECEL

  def get_lead(self, sm):

    radar = sm['radarState']
    model = sm['model']

    if radar.leadOne.status and radar.leadOne.modelProb >= 0.5:
      return radar.leadOne

    try:
      radar = log.RadarState.LeadData.new_message()
      radar.leadOne.status = 1
      radar.leadOne.modelProb = model.lead.prob
      radar.leadOne.dRel = model.lead.dist - RADAR_TO_CAMERA
      radar.leadOne.vRel = model.lead.relVel
      radar.leadOne.yRel = model.lead.relY

      #radar.leadTwo.status = 1
      #radar.leadTwo.modelProb = model.leadFuture.prob
      #radar.leadTwo.dRel = model.leadFuture.dist - RADAR_TO_CAMERA
      #radar.leadTwo.vRel = model.leadFuture.relVel
      #radar.leadTwo.yRel = model.leadFuture.relY

      return radar.leadOne
    except:
      pass

    return None

  def cal_acc(self, apply_accel, CS, clu11_speed, sm):

    cruise_gap = clip(CS.cruise_gap, 1., 4.)

    override_acc = 0.
    v_ego = clu11_speed * CV.KPH_TO_MS
    op_accel = apply_accel * CV.MS_TO_KPH

    lead = self.get_lead(sm)
    if lead is None or lead.modelProb < 0.5:
      accel = op_accel
    else:

      d = lead.dRel - 5.

      if 0. < d < -lead.vRel * (7.7 + cruise_gap) * 2. and lead.vRel < -1.:  #여기서 부터 254라인까지 선행차가 있을 때 감속 로직, 내용이 있으나 고치면 생각보다 처음부터 해야할게 많아서 안건드리는게 낫다
        t = d / lead.vRel * 0.98
        acc = -(lead.vRel / t) * CV.MS_TO_KPH * 1.8
        override_acc = acc
        accel = (op_accel + acc) / 2.
      else:                #여기서 부터 선행차가 있을 때 가속 로직.        
        accel = op_accel * interp(clu11_speed, [0., 25., 50., 51., 60., 100.], [2.4, 2.9, 1.7, 1.65, 1.4, 1.0])    #보간법 0일 때 엑셀값(가속도) 2.4배 나머지는 각 구간의 보간법
#        if 35 > lead.dRel > 15:   257부터 261 오류 작동안함.
#          if clu11_speed < 50:
#            accel = op_accel * 2.5               
#        else:
#          accel = op_accel * interp(clu11_speed, [50., 60., 100.], [1.75, 1.4, 1.0])
    if accel > 0.:    # 전체 가속로직 선행차가 있을 땐 256라인 가속을 받고 여기서 한번더 가속을 받는다. 선행차가 없을 때는 여기서만 가속을 받는다. 자연스러운 가속은 선행차가 있을 때 가속을 세게해야한다. 
      accel *= self.accel_gain * interp(clu11_speed, [35., 60., 100.], [1.5, 1.25, 1.215])
    else:
      accel *= self.decel_gain * interp(clu11_speed, [70., 75.], [1.79287, 1.8])    #전체 감속로직 보간70킬로까지 감속을 너무 세게 하면 자연스러운 거리 유지가 힘듬.

    return clip(accel, -LIMIT_DECEL, LIMIT_ACCEL), override_acc

  def get_curve_speed(self, sm, v_ego):

    if not self.slow_on_curves:
      return 255.

    if len(sm['model'].path.poly):
      path = list(sm['model'].path.poly)

      # Curvature of polynomial https://en.wikipedia.org/wiki/Curvature#Curvature_of_the_graph_of_a_function
      # y = a x^3 + b x^2 + c x + d, y' = 3 a x^2 + 2 b x + c, y'' = 6 a x + 2 b
      # k = y'' / (1 + y'^2)^1.5
      # TODO: compute max speed without using a list of points and without numpy
      y_p = 3 * path[0] * self.path_x**2 + 2 * path[1] * self.path_x + path[2]
      y_pp = 6 * path[0] * self.path_x + 2 * path[1]
      curv = y_pp / (1. + y_p**2)**1.5

      a_y_max = 2.975 - v_ego * 0.0375  # ~1.85 @ 75mph, ~2.6 @ 25mph
      v_curvature = np.sqrt(a_y_max / np.clip(np.abs(curv), 1e-4, None))
      model_speed = np.min(v_curvature) * self.curvature_gain
      model_speed = max(32. * CV.KPH_TO_MS, model_speed) # Don't slow down below 32km/h
    else:
      model_speed = 255.

    return model_speed

  @staticmethod
  def update_cruise_buttons(controls, CS):

    car_set_speed = CS.cruiseState.speed * CV.MS_TO_KPH
    is_cruise_enabled = car_set_speed != 0 and car_set_speed != 255 and CS.cruiseState.enabled and controls.CP.enableCruise

    if is_cruise_enabled:
      if controls.CC.sccSmoother.state == CruiseState.STOCK:
        controls.v_cruise_kph = CS.cruiseState.speed * CV.MS_TO_KPH
      else:
        controls.v_cruise_kph = SccSmoother.update_v_cruise(controls.v_cruise_kph, CS.buttonEvents, controls.enabled, controls.is_metric)
    else:
      controls.v_cruise_kph = 0

    if controls.is_cruise_enabled != is_cruise_enabled:
      controls.is_cruise_enabled = is_cruise_enabled

      if controls.is_cruise_enabled:
        controls.v_cruise_kph = CS.cruiseState.speed * CV.MS_TO_KPH
      else:
        controls.v_cruise_kph = 0

      controls.LoC.reset(v_pid=CS.vEgo)

    #controls.v_cruise_kph = road_speed_limiter_get_max_speed(CS, controls.v_cruise_kph)
    controls.cruiseOpMaxSpeed = controls.v_cruise_kph

  @staticmethod
  def update_v_cruise(v_cruise_kph, buttonEvents, enabled, metric):

    global ButtonCnt, LongPressed, ButtonPrev
    if enabled:
      if ButtonCnt:
        ButtonCnt += 1
      for b in buttonEvents:
        if b.pressed and not ButtonCnt and (b.type == ButtonType.accelCruise or b.type == ButtonType.decelCruise):
          ButtonCnt = 1
          ButtonPrev = b.type
        elif not b.pressed and ButtonCnt:
          if not LongPressed and b.type == ButtonType.accelCruise:
            v_cruise_kph += 1 if metric else 1 * CV.MPH_TO_KPH
          elif not LongPressed and b.type == ButtonType.decelCruise:
            v_cruise_kph -= 1 if metric else 1 * CV.MPH_TO_KPH
          LongPressed = False
          ButtonCnt = 0
      if ButtonCnt > 70:
        LongPressed = True
        V_CRUISE_DELTA = V_CRUISE_DELTA_KM if metric else V_CRUISE_DELTA_MI
        if ButtonPrev == ButtonType.accelCruise:
          v_cruise_kph += V_CRUISE_DELTA - v_cruise_kph % V_CRUISE_DELTA
        elif ButtonPrev == ButtonType.decelCruise:
          v_cruise_kph -= V_CRUISE_DELTA - -v_cruise_kph % V_CRUISE_DELTA
        ButtonCnt %= 70
      v_cruise_kph = clip(v_cruise_kph, MIN_SET_SPEED, MAX_SET_SPEED)

    return v_cruise_kph



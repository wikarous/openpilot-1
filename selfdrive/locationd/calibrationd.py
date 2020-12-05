#!/usr/bin/env python3
'''
This process finds calibration values. More info on what these calibration values
are can be found here https://github.com/commaai/openpilot/tree/master/common/transformations
While the roll calibration is a real value that can be estimated, here we assume it's zero,
and the image input into the neural network is not corrected for roll.
'''

import os
import copy
import json
import numpy as np
import cereal.messaging as messaging
from selfdrive.config import Conversions as CV
from selfdrive.swaglog import cloudlog
from common.params import Params, put_nonblocking
from common.transformations.model import model_height
from common.transformations.camera import get_view_frame_from_road_frame
from common.transformations.orientation import rot_from_euler, euler_from_rot

MIN_SPEED_FILTER = 15 * CV.MPH_TO_MS
MAX_VEL_ANGLE_STD = np.radians(0.25)
MAX_YAW_RATE_FILTER = np.radians(2)  # per second

# This is at model frequency, blocks needed for efficiency
SMOOTH_CYCLES = 400
BLOCK_SIZE = 100
INPUTS_NEEDED = 5   # Minimum blocks needed for valid calibration
INPUTS_WANTED = 50   # We want a little bit more than we need for stability
MAX_ALLOWED_SPREAD = np.radians(2)
RPY_INIT = np.array([0.0,0.0,0.0])

# These values are needed to accomodate biggest modelframe
PITCH_LIMITS = np.array([-0.09074112085129739, 0.14907572052989657])
YAW_LIMITS = np.array([-0.06912048084718224, 0.06912048084718235])
DEBUG = os.getenv("DEBUG") is not None


class Calibration:
  UNCALIBRATED = 0
  CALIBRATED = 1
  INVALID = 2


def is_calibration_valid(rpy):
  return (PITCH_LIMITS[0] < rpy[1] < PITCH_LIMITS[1]) and (YAW_LIMITS[0] < rpy[2] < YAW_LIMITS[1])


def sanity_clip(rpy):
  if np.isnan(rpy).any():
    rpy = RPY_INIT
  return np.array([rpy[0],
                   np.clip(rpy[1], PITCH_LIMITS[0] - .005, PITCH_LIMITS[1] + .005),
                   np.clip(rpy[2], YAW_LIMITS[0] - .005, YAW_LIMITS[1] + .005)])


class Calibrator():
  def __init__(self, param_put=False):
    self.param_put = param_put

    # Read saved calibration
    calibration_params = Params().get("CalibrationParams")

    rpy_init = RPY_INIT
    valid_blocks = 0

    if param_put and calibration_params:
      try:
        calibration_params = json.loads(calibration_params)
        rpy_init = calibration_params["calib_radians"]
        valid_blocks = calibration_params['valid_blocks']
      except Exception:
        cloudlog.exception("CalibrationParams file found but error encountered")
    self.reset(rpy_init, valid_blocks)
    self.update_status()

  def reset(self, rpy_init=RPY_INIT, valid_blocks=0, smooth_from=None):
    if not np.isfinite(rpy_init).all():
        self.rpy = copy.copy(RPY_INIT)
    else:
      self.rpy = rpy_init
    if not np.isfinite(valid_blocks) or valid_blocks < 0:
        self.valid_blocks = 0
    else:
      self.valid_blocks = valid_blocks
    self.rpys = np.tile(self.rpy, (INPUTS_WANTED, 1))

    self.idx = 0
    self.block_idx = 0
    self.v_ego = 0

    if smooth_from is None:
      self.old_rpy = RPY_INIT
      self.old_rpy_weight = 0.0
    else:
      self.old_rpy = smooth_from
      self.old_rpy_weight = 1.0

  def update_status(self):
    if self.valid_blocks > 0:
      max_rpy_calib = np.array(np.max(self.rpys[:self.valid_blocks], axis=0))
      min_rpy_calib = np.array(np.min(self.rpys[:self.valid_blocks], axis=0))
      self.calib_spread = np.abs(max_rpy_calib - min_rpy_calib)
    else:
      self.calib_spread = np.zeros(3)

    if self.valid_blocks < INPUTS_NEEDED:
      self.cal_status = Calibration.UNCALIBRATED
    elif is_calibration_valid(self.rpy):
      self.cal_status = Calibration.CALIBRATED
    else:
      self.cal_status = Calibration.INVALID

    # If spread is too high, assume mounting was changed and reset to last block.
    # Make the transition smooth. Abrupt transistion are not good foor feedback loop through supercombo model.
    if max(self.calib_spread) > MAX_ALLOWED_SPREAD and self.cal_status == Calibration.CALIBRATED:
      self.reset(self.rpys[self.block_idx - 1], valid_blocks=INPUTS_NEEDED, smooth_from=self.rpy)

    write_this_cycle = (self.idx == 0) and (self.block_idx % (INPUTS_WANTED//5) == 5)
    if self.param_put and write_this_cycle:
      # TODO: this should use the liveCalibration struct from cereal
      cal_params = {"calib_radians": list(self.rpy),
                    "valid_blocks": int(self.valid_blocks)}
      put_nonblocking("CalibrationParams", json.dumps(cal_params).encode('utf8'))

  def handle_v_ego(self, v_ego):
    self.v_ego = v_ego

  def get_smooth_rpy(self):
    if self.old_rpy_weight > 0:
      return self.old_rpy_weight * self.old_rpy + (1.0 - self.old_rpy_weight) * self.rpy
    else:
      return self.rpy

  def handle_cam_odom(self, trans, rot, trans_std, rot_std):
    self.old_rpy_weight = min(0.0, self.old_rpy_weight - 1/SMOOTH_CYCLES)

    straight_and_fast = ((self.v_ego > MIN_SPEED_FILTER) and (trans[0] > MIN_SPEED_FILTER) and (abs(rot[2]) < MAX_YAW_RATE_FILTER))
    certain_if_calib = ((np.arctan2(trans_std[1], trans[0]) < MAX_VEL_ANGLE_STD) or
                        (self.valid_blocks < INPUTS_NEEDED))
    if straight_and_fast and certain_if_calib:
      observed_rpy = np.array([0,
                               -np.arctan2(trans[2], trans[0]),
                               np.arctan2(trans[1], trans[0])])
      new_rpy = euler_from_rot(rot_from_euler(self.get_smooth_rpy()).dot(rot_from_euler(observed_rpy)))
      new_rpy = sanity_clip(new_rpy)

      self.rpys[self.block_idx] = (self.idx*self.rpys[self.block_idx] + (BLOCK_SIZE - self.idx) * new_rpy) / float(BLOCK_SIZE)
      self.idx = (self.idx + 1) % BLOCK_SIZE
      if self.idx == 0:
        self.block_idx += 1
        self.valid_blocks = max(self.block_idx, self.valid_blocks)
        self.block_idx = self.block_idx % INPUTS_WANTED
      if self.valid_blocks > 0:
        self.rpy = np.mean(self.rpys[:self.valid_blocks], axis=0)


      self.update_status()

      return new_rpy
    else:
      return None

  def send_data(self, pm):
    smooth_rpy = self.get_smooth_rpy()
    extrinsic_matrix = get_view_frame_from_road_frame(0, smooth_rpy[1], smooth_rpy[2], model_height)

    cal_send = messaging.new_message('liveCalibration')
    cal_send.liveCalibration.validBlocks = self.valid_blocks
    cal_send.liveCalibration.calStatus = self.cal_status
    cal_send.liveCalibration.calPerc = min(100 * (self.valid_blocks * BLOCK_SIZE + self.idx) // (INPUTS_NEEDED * BLOCK_SIZE), 100)
    cal_send.liveCalibration.extrinsicMatrix = [float(x) for x in extrinsic_matrix.flatten()]
    cal_send.liveCalibration.rpyCalib = [float(x) for x in smooth_rpy]
    cal_send.liveCalibration.rpyCalibSpread = [float(x) for x in self.calib_spread]

    pm.send('liveCalibration', cal_send)


def calibrationd_thread(sm=None, pm=None):
  if sm is None:
    sm = messaging.SubMaster(['cameraOdometry', 'carState'], poll=['cameraOdometry'])

  if pm is None:
    pm = messaging.PubMaster(['liveCalibration'])

  calibrator = Calibrator(param_put=True)

  while 1:
    timeout = 0 if sm.frame == -1 else 100
    sm.update(timeout)

    if sm.updated['cameraOdometry']:
      calibrator.handle_v_ego(sm['carState'].vEgo)
      new_rpy = calibrator.handle_cam_odom(sm['cameraOdometry'].trans,
                                           sm['cameraOdometry'].rot,
                                           sm['cameraOdometry'].transStd,
                                           sm['cameraOdometry'].rotStd)

      if DEBUG and new_rpy is not None:
        print('got new rpy', new_rpy)

    # 4Hz driven by cameraOdometry
    if sm.frame % 5 == 0:
      calibrator.send_data(pm)


def main(sm=None, pm=None):
  calibrationd_thread(sm, pm)


if __name__ == "__main__":
  main()

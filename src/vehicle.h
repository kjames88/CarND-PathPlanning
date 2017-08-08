
#pragma once

#include <deque>

class Vehicle {
 public:
  constexpr static const double timestep_ = 0.02;
 Vehicle() : id_(-1), lane_(-1), s_(0.0), d_(0.0), velocity_(0.0), acceleration_(0.0) {};
 Vehicle(int id, int lane, double s, double d, double v) : id_(id), lane_(lane), s_(s), d_(d), velocity_(v), acceleration_(0.0) {};
  void set_lane(int lane) {lane_ = lane;}
  void set_s(double s) {s_ = s;}
  void set_d(double d) {d_ = d;}
  void set_velocity(double v) {
    velocity_ = v;
    push_v(v);
  }
  void set_acceleration(double a) {acceleration_ = a;}
  int get_lane() {return lane_;}
  bool uses_lane(int lane) {
    if (lane == lane_) {
      return true;
    }
    // check for lane change in progress
    // 1m off center takes the next lane as well as current
    double center = d_;
    while (center >= 4.0) {
      center -= 4.0;
    }
    if (lane == (lane_ + 1)) {
      if (center >= 3.0) {
        return true;
      }
    } else if (lane == (lane_ - 1)) {
      if (center <= 1.0) {
        return true;
      }
    }
    return false;
  }
  double get_s() {return s_;}
  double get_velocity() {return velocity_;}
  double get_acceleration() {return acceleration_;}
 private:
  int id_;
  int lane_;
  double s_;
  double d_;
  double velocity_;
  std::deque<double> v_samples_;
  double acceleration_;

  void push_v(double v) {
    v_samples_.push_back(v);
    if (v_samples_.size() > 10) {
      v_samples_.pop_front();
    }
    acceleration_ = (v_samples_.back() - v_samples_.front()) / ((double) v_samples_.size() * timestep_);
  }
};

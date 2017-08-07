#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <map>
#include <time.h>
#include <assert.h>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "Eigen-3.3/Eigen/LU"
#include "json.hpp"
#include "vehicle.h"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// The max s value before wrapping around the track back to 0
double max_s = 6945.554;

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
  return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

  double closestLen = 100000; //large number
  int closestWaypoint = 0;

  for(int i = 0; i < maps_x.size(); i++)
    {
      double map_x = maps_x[i];
      double map_y = maps_y[i];
      double dist = distance(x,y,map_x,map_y);
      if(dist < closestLen)
        {
          closestLen = dist;
          closestWaypoint = i;
        }

    }

  return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

  int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

  double map_x = maps_x[closestWaypoint];
  double map_y = maps_y[closestWaypoint];

  double heading = atan2( (map_y-y),(map_x-x) );

  double angle = abs(theta-heading);

  if(angle > pi()/4)
    {
      closestWaypoint++;
    }

  return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
  int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

  int prev_wp;
  prev_wp = next_wp-1;
  if(next_wp == 0)
    {
      prev_wp  = maps_x.size()-1;
    }

  double n_x = maps_x[next_wp]-maps_x[prev_wp];
  double n_y = maps_y[next_wp]-maps_y[prev_wp];
  double x_x = x - maps_x[prev_wp];
  double x_y = y - maps_y[prev_wp];

  // find the projection of x onto n
  double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
  double proj_x = proj_norm*n_x;
  double proj_y = proj_norm*n_y;

  double frenet_d = distance(x_x,x_y,proj_x,proj_y);

  //see if d value is positive or negative by comparing it to a center point

  double center_x = 1000-maps_x[prev_wp];
  double center_y = 2000-maps_y[prev_wp];
  double centerToPos = distance(center_x,center_y,x_x,x_y);
  double centerToRef = distance(center_x,center_y,proj_x,proj_y);

  if(centerToPos <= centerToRef)
    {
      frenet_d *= -1;
    }

  // calculate s value
  double frenet_s = 0;
  for(int i = 0; i < prev_wp; i++)
    {
      frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
    }

  frenet_s += distance(0,0,proj_x,proj_y);

  return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
  int prev_wp = -1;

  while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
    {
      prev_wp++;
    }

  int wp2 = (prev_wp+1)%maps_x.size();

  double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
  // the x,y,s along the segment
  double seg_s = (s-maps_s[prev_wp]);

  double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
  double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

  double perp_heading = heading-pi()/2;

  double x = seg_x + d*cos(perp_heading);
  double y = seg_y + d*sin(perp_heading);

  return {x,y};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
//   - extended withd dx,dy since the initial method drifted in d using spline fit hires map
vector<double> getXY(double s, double d,
                     vector<double>& maps_s, vector<double>& maps_x, vector<double>& maps_y,
                     vector<double>& maps_dx, vector<double>& maps_dy) {
  
  assert(s >= 0.0 && s <= max_s);
  int prev_wp = -1;

  while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
    {
      prev_wp++;
    }
  int wp2 = (prev_wp+1)%maps_x.size();
  double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
  double delta = s - maps_s[prev_wp];
  
  double x = maps_x[prev_wp] + (delta * cos(heading)) + (d * maps_dx[prev_wp]);
  double y = maps_y[prev_wp] + (delta * sin(heading)) + (d * maps_dy[prev_wp]);
  return {x,y};
}

vector<double> getXY(double s, double d,
                     tk::spline& sx, tk::spline& sy, tk::spline& sdx, tk::spline& sdy) {
  if (s < 0.0 || s > max_s) {
    std::cout << "s=" << s << " out of range!" << std::endl;
  }
  assert(s >= 0.0 && s <= max_s);
  double x = sx(s) + (d * sdx(s));
  double y = sy(s) + (d * sdy(s));
  return {x,y};
}

double get_map_s(double s) {
  double s_map = s;
  while (s_map > max_s) {
    s_map -= max_s;
  }
  return s_map;
}

std::vector<double> solve_quintic(std::vector<double>& poly, double t) {
  //std::cout << "solve_quintic t=" << t << " poly " << poly[0] << "," << poly[1] << "," << poly[2] << "," << poly[3] << "," << poly[4] << "," << poly[5] << std::endl;
  double s = poly[0] + (poly[1] * t) + (poly[2] * pow(t,2)) + (poly[3] * pow(t,3)) + (poly[4] * pow(t,4)) + (poly[5] * pow(t,5));
  double v = poly[1] + (2.0 * poly[2] * t) + (3.0 * poly[3] * pow(t,2)) + (4.0 * poly[4] * pow(t,3)) + (5.0 * poly[5] * pow(t,4));
  double a = (2.0 * poly[2]) + (6.0 * poly[3] * t) + (12.0 * poly[4] * pow(t,2)) + (20.0 * poly[5] * pow(t,3));
  return {s,v,a};
}

double mph_to_mps(double mph) {
  double mps = mph / 3600.0;  // miles/sec
  mps *= (0.0254 * 12.0 * 5280.0);  // m/mile
  return mps;
}

// check the next 50 timesteps for the path for speed, acceleration, and jerk
bool check_path(std::vector<double> s_poly, std::vector<double> d_poly,
                tk::spline& sx, tk::spline& sy, tk::spline& sdx, tk::spline& sdy) {
  bool fail = false;
  double max_speed = mph_to_mps(50.0);
  double max_acc = 10.0;
  double max_jerk = 10.0;
  double x_prev[3] = {0.0, 0.0, 0.0};
  double y_prev[3] = {0.0, 0.0, 0.0};
  double v_prev = 0.0;
  double a_prev = 0.0;
  for (int i=0; i<50; i++) {
    double t = (double) i * 0.02;
    auto s_soln = solve_quintic(s_poly, t);
    auto d_soln = solve_quintic(d_poly, t);
    double s = get_map_s(s_soln[0]);
    double d = d_soln[0];
    auto dbl_vec = getXY(s, d, sx, sy, sdx, sdy);
    double x = dbl_vec[0];
    double y = dbl_vec[1];
    if (i>0) {
      // check speed
      double v = sqrt(pow(x-x_prev[0],2) + pow(y-y_prev[0],2)) / 0.02;
      if (v > max_speed) {
        std::cout << "FAIL v=" << v << " idx " << i << std::endl;
        fail = true;
        return false;
      }
      if (i>1) {
        // check acceleration
        double a = (v-v_prev) / 0.02;
        if (a > max_acc) {
          std::cout << "FAIL a=" << a << std::endl;
          fail = true;
          return false;
        }
        if (i>2) {
          // check jerk
          double j = (a - a_prev) / 0.02;
          if (j > max_jerk) {
            // std::cout << "FAIL j=" << j << std::endl;
            // fail = true;
            // return false;
          }
        }
        a_prev = a;
      }
      v_prev = v;
    }
    for (int j=2; j>0; j--) {
      x_prev[j] = x_prev[j-1];
      y_prev[j] = y_prev[j-1];
    }
    x_prev[0] = x;
    y_prev[0] = y;
  }
  return !fail;
}

// check the next 50 timesteps for the path for speed, acceleration, and jerk
bool check_path(std::vector<double> next_x, std::vector<double> next_y) {
  bool fail = false;
  assert(next_x.size() > 0);
  int steps = next_x.size() > 50 ? 50 : next_x.size();
  double max_speed = mph_to_mps(50.0);
  double max_acc = 10.0;
  double max_jerk = 10.0;
  double x_prev[3] = {0.0, 0.0, 0.0};
  double y_prev[3] = {0.0, 0.0, 0.0};
  double v_prev = 0.0;
  double a_prev = 0.0;
  for (int i=0; i<50; i++) {
    double t = (double) i * 0.02;
    double x = next_x[i];
    double y = next_y[i];
    if (i>0) {
      // check speed
      double v = sqrt(pow(x-x_prev[0],2) + pow(y-y_prev[0],2)) / 0.02;
      if (v > max_speed) {
        std::cout << "FAIL v=" << v << " idx " << i << std::endl;
        fail = true;
        return false;
      }
      if (i>1) {
        // check acceleration
        double a = (v-v_prev) / 0.02;
        if (a > max_acc) {
          std::cout << "FAIL a=" << a << std::endl;
          fail = true;
          return false;
        }
        if (i>2) {
          // check jerk
          double j = (a - a_prev) / 0.02;
          if (j > max_jerk) {
            // std::cout << "FAIL j=" << j << ", a=" << a << " a_prev=" << a_prev << std::endl;
            // fail = true;
            // return false;
          }
        }
        a_prev = a;
      }
      v_prev = v;
    }
    for (int j=2; j>0; j--) {
      x_prev[j] = x_prev[j-1];
      y_prev[j] = y_prev[j-1];
    }
    x_prev[0] = x;
    y_prev[0] = y;
  }
  return !fail;
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
    istringstream iss(line);
    double x;
    double y;
    float s;
    float d_x;
    float d_y;
    iss >> x;
    iss >> y;
    iss >> s;
    iss >> d_x;
    iss >> d_y;
    map_waypoints_x.push_back(x);
    map_waypoints_y.push_back(y);
    map_waypoints_s.push_back(s);
    map_waypoints_dx.push_back(d_x);
    map_waypoints_dy.push_back(d_y);
  }

  // Create a wrap-around waypoint for the splines to follow from the track end to the track start
  //   x,y,dx,dy are all those of the first waypoint on the track
  //   for s we'll create an out of range s coordinate at the first waypoint
  double s_first = map_waypoints_s.front();
  double s_last = map_waypoints_s.back();
  double s_extended = s_first + max_s;
  map_waypoints_s.push_back(s_extended);
  map_waypoints_x.push_back(map_waypoints_x.front());
  map_waypoints_y.push_back(map_waypoints_y.front());
  map_waypoints_dx.push_back(map_waypoints_dx.front());
  map_waypoints_dy.push_back(map_waypoints_dy.front());
  
  // following what others are mentioning on forum and slack, instead of spline x,y waypoints,
  //   spline s,x and s,y and then mapping from s should be simple
  tk::spline sx;
  tk::spline sy;
  tk::spline sdx;
  tk::spline sdy;
  sx.set_points(map_waypoints_s, map_waypoints_x);
  sy.set_points(map_waypoints_s, map_waypoints_y);
  sdx.set_points(map_waypoints_s, map_waypoints_dx);
  sdy.set_points(map_waypoints_s, map_waypoints_dy);
  
  std::map<int, Vehicle> vehicles;
  std::vector<double> prev_x_vals;
  std::vector<double> prev_s_vals;
  std::vector<double> prev_d_vals;
  int car_lane = 0;
  bool lane_change_state = false;
  int lane_change_lane = 0;
  double lane_change_speed = 0.0;
  int msg_cnt = 0;
  double max_ratio = 1.0;      // x,y distance to s distance ratio for regulating s_dot
  h.onMessage([&max_ratio, &vehicles, &car_lane,
               &lane_change_state, &lane_change_lane, &lane_change_speed,
               &msg_cnt, &prev_x_vals, &prev_s_vals, &prev_d_vals,
               &sx, &sy, &sdx, &sdy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                                     uWS::OpCode opCode) {
                // "42" at the start of the message means there's a websocket message event.
                // The 4 signifies a websocket message
                // The 2 signifies a websocket event
                //auto sdata = string(data).substr(0, length);
                if (length && length > 2 && data[0] == '4' && data[1] == '2') {

                  auto s = hasData(data);

                  if (s != "") {
                    auto j = json::parse(s);
        
                    string event = j[0].get<string>();
        
                    if (event == "telemetry") {
                      // j[1] is the data JSON object
          
                      // Main car's localization Data
                      double car_x = j[1]["x"];
                      double car_y = j[1]["y"];
                      double car_s = j[1]["s"];
                      double car_d = j[1]["d"];
                      double car_yaw = j[1]["yaw"];
                      double car_speed = j[1]["speed"];

                      // Previous path data given to the Planner
                      auto previous_path_x = j[1]["previous_path_x"];
                      auto previous_path_y = j[1]["previous_path_y"];
                      // Previous path's end s and d values 
                      double end_path_s = j[1]["end_path_s"];
                      double end_path_d = j[1]["end_path_d"];

                      // Sensor Fusion Data, a list of all other cars on the same side of the road.
                      //  Format:  [id,x,y,vx,vy,s,d]
                      auto sensor_fusion = j[1]["sensor_fusion"];

                      //std::cout << "sensor_fusion: " << sensor_fusion << std::endl;
                      for (int i=0; i<sensor_fusion.size(); i++) {
                        int id = sensor_fusion[i][0];
                        double obstacle_d = sensor_fusion[i][6];
                        int lane = obstacle_d / 4.0;
                        double obstacle_s = sensor_fusion[i][5];
                        double obstacle_vx = sensor_fusion[i][3];
                        double obstacle_vy = sensor_fusion[i][4];
                        double obstacle_v = sqrt(pow(obstacle_vx,2) + pow(obstacle_vy,2));
                        auto it = vehicles.find(id);
                        if (it == vehicles.end()) {
                          Vehicle v(id, lane, obstacle_s, obstacle_v);
                          it = vehicles.insert(std::pair<int, Vehicle> (id, v)).first;
                        } else {
                          it->second.set_lane(lane);
                          it->second.set_s(obstacle_s);
                          it->second.set_velocity(obstacle_v);
                        }
                      }

                      // check for a car in front of me to follow
                      //   find any vehicle in my lane (assume obstacle vehicles don't drive on the lines)
                      //   sort by closest distance
                      double T = 3.0;
                      assert(max_ratio >= 1.0);
                      double maximum_speed = mph_to_mps(48.5 / max_ratio);
                      double car_v = mph_to_mps(car_speed);
                      double car_a = 0.0;
                      double car_vd = 0.0;
                      double car_ad = 0.0;
                      bool have_target_vehicle = false;
                      int target_vehicle = -1;
                      double min_dist = 1e6;
                      double min_velocity = 0.0;
                      int min_id = -1;
                      
                      // Either complete a lane change or look for someone to follow
                      if (lane_change_state == false) {
                        for (auto it=vehicles.begin(); it!=vehicles.end(); it++) {
                          if (it->second.get_lane() == car_lane) {
                            double dist = it->second.get_s() - car_s;
                            if (dist >= 0 && dist < min_dist) {
                              min_dist = dist;
                              min_id = it->first;
                              min_velocity = it->second.get_velocity();
                            }
                          }
                        }
                        if (min_id >= 0 && min_dist < 50.0) {
                          auto it = vehicles.find(min_id);
                          assert (it != vehicles.end());
                          // std::cout << "LEADER id " << min_id << " distance " << min_dist << "m velocity "
                          //           << it->second.get_velocity() << "m/s" << std::endl;
                          target_vehicle = min_id;
                          
                          // FIXME unify duplicated code
                          
                          double s_lv = vehicles[target_vehicle].get_s();             // lead vehicle initial s
                          double v_lv = vehicles[target_vehicle].get_velocity();      // lead vehicle velocity
                          double a_lv = vehicles[target_vehicle].get_acceleration();  // lead vehicle acceleration
                          double sf_lv = s_lv + (v_lv * T) + (0.5 * a_lv * pow(T,2)); // may be out of range
                          double sf_max = get_map_s(sf_lv - 15.0);  // don't get closer than 15m (could be made speed dependent)
                          if (get_map_s(car_s + (car_v * T) + (0.5 * car_a * pow(T,2))) >= sf_max) {
                            have_target_vehicle = true;
                          }
                        }
                      } else {
                        // keep the parameters established when we decided to change lanes
                      }
                      if (have_target_vehicle) {
                        // following someone slower:  see if we can drive faster by passing in a neighboring lane
                        bool try_lane[3] = {car_lane == 1, car_lane == 2 || car_lane == 0, car_lane == 1};
                        bool block_lane[3] = {false, false, false};
                        double lane_vmax[3] = {maximum_speed, maximum_speed, maximum_speed};
                        for (auto it=vehicles.begin(); it!=vehicles.end(); it++) {
                          int veh_lane = it->second.get_lane();
                          if (veh_lane >= 0 && veh_lane < 3 && try_lane[veh_lane]) {
                            double s0 = it->second.get_s();
                            double s1 = get_map_s(s0 + (T * it->second.get_velocity())
                                                  + (0.5 * it->second.get_acceleration() * pow(T,2)));
                            // block the lane if a car is within buffer distance in front or behind
                            double dist = abs(car_s - s0);
                            bool ahead = (car_s > s0);
                            double car_s1 = get_map_s(car_s + (car_v * T) + (0.5 * car_a * pow(T,2)));
                            if ((ahead && (dist < 5.0)) || (!ahead && (dist < 10.0))) {
                              block_lane[veh_lane] = true;
                              // std::cout << "veh " << it->first << " s=" << s0 << " car_s=" << car_s <<
                              //   " blocks lane " << veh_lane << " (1)" << std::endl;
                            } else if (!ahead && ((s1 - car_s1) < 20.0)) {
                              block_lane[veh_lane] = true;
                              // std::cout << "veh " << it->first << " s=" << s0 << "," << s1 << " car_s1=" << car_s1 <<
                              //   " blocks lane " << veh_lane << " (1.5)" << std::endl;
                            } else {
                              // block the lane if a car will cross our position at velocity
                              //   - also ignore a lane with a slower car in front
                              if (s0 < car_s) {
                                if (s1 > get_map_s(car_s + (car_v * T) + (0.5 * car_a * pow(T,2)))) {
                                  block_lane[veh_lane] = true;
                                  // std::cout << "veh " << it->first << " blocks lane " << veh_lane << " (2)" << std::endl;
                                }
                              } else {
                                if (dist < 50.0 && car_v > it->second.get_velocity()) {
                                  block_lane[veh_lane] = true;
                                  // std::cout << "veh " << it->first << " blocks lane " << veh_lane << " (3)" << std::endl;
                                } else {
                                  if (dist < 50.0) {
                                    lane_vmax[veh_lane] = it->second.get_velocity();
                                  }
                                }
                              }
                            }
                          }
                        }
                        int max_lane = -1;
                        double max_v = 0.0;
                        for (int i=0; i<3; i++) {
                          if (try_lane[i] && !block_lane[i]) {
                            max_lane = i;
                            max_v = lane_vmax[i];
                          }
                        }
                        if (max_lane >= 0) {
                          // std::cout << "we can move into lane " << max_lane << " v=" << max_v << std::endl;
                          lane_change_state = true;
                          lane_change_lane = max_lane;
                          lane_change_speed = max_v;
                          double max_lane_change_speed = 0.9 * (mph_to_mps(50.0));  // account for curvature and d component
                          if (lane_change_speed > max_lane_change_speed) {
                            lane_change_speed = max_lane_change_speed;
                          }
                        }
                      }
                      
                      json msgJson;

                      vector<double> next_x_vals;
                      vector<double> next_y_vals;
                      vector<double> next_s_vals;
                      vector<double> next_d_vals;
                      int copy_path_cnt = 5;
                      if (car_v > 0) {

                        // determine the number of points popped
                        //   then we can use prior information for the start of the next segment
                        int popped = prev_x_vals.size() - previous_path_x.size();
                        assert(prev_x_vals.size() == prev_s_vals.size());
                        assert(prev_x_vals.size() == prev_d_vals.size());
                        assert(popped > 0);
                      
                        // stitch in the beginning of the previous path and then project from there for continuity
                        double init_x = previous_path_x[0];
                        double init_y = previous_path_y[0];
                        int prev=0;
                        //while (distance(init_x, init_y, previous_path_x[prev], previous_path_y[prev]) < 2.0
                        while (prev < copy_path_cnt
                               && prev < previous_path_x.size() && previous_path_x[prev] != NULL) {
                          next_x_vals.push_back(previous_path_x[prev]);
                          next_y_vals.push_back(previous_path_y[prev]);
                          prev++;
                        }
                        copy_path_cnt = prev;
                        // keep copied s,d in prev vectors; sync'd with next_x_vals and next_y_vals
                        for (int i=0; i<copy_path_cnt; i++) {
                          next_s_vals.push_back(prev_s_vals[popped+i]);
                          next_d_vals.push_back(prev_d_vals[popped+i]);
                        }
                        assert(prev >= 3);

                        if (prev > 0) {
                          // push only as many next points as we want to stitch into the new path so the s,d values can be used
                          //   - also keep the final v,a from the previous polynomial
                          int idx_back = next_s_vals.size() - 1;
                          car_s = next_s_vals[idx_back];
                          car_v = (next_s_vals[idx_back] - next_s_vals[idx_back-1]) / 0.02;
                          car_a = (car_v - ((next_s_vals[idx_back-1] - next_s_vals[idx_back-2]) / 0.02)) / 0.02;

                          car_d = next_d_vals[idx_back];
                          car_vd = (next_d_vals[idx_back] - next_d_vals[idx_back-1]) / 0.02;
                          car_ad = (car_vd - ((next_d_vals[idx_back-1] - next_d_vals[idx_back-2]) / 0.02)) / 0.02;
                        }                        
                        // std::cout << "popped=" << popped << ", copy_path_cnt=" << copy_path_cnt << " plot from car_s=" << car_s
                        //           << " car_v=" << car_v << " car_a=" << car_a << std::endl;                        
                        
                      } else {

                        std::cout << "Car is stopped" << std::endl;
                        
                        copy_path_cnt = 0;
                        car_lane = car_d / 4.0;  // start out in the lane set by the simulator init
                      }

                      double si = car_s;
                      double si_dot = car_v;
                      double si_dot_dot = car_a;
                      double di = car_d;
                      double di_dot = car_vd;
                      double di_dot_dot = car_ad;

                      std::vector<double> s_poly;
                      std::vector<double> d_poly;
                      vector<double> x_vals_raw;
                      vector<double> y_vals_raw;
                      vector<double> s_vals_raw;
                      vector<double> d_vals_raw;

                      bool searching = true;
                      int search_cnt = 0;
                      while (searching) {
                      
                        // Using the quintic polynomial solver from Trajectory Generation
                        //   Or a quadric alternative to get the car moving suggested in the Werling paper

                        x_vals_raw.clear();
                        y_vals_raw.clear();
                        s_vals_raw.clear();
                        d_vals_raw.clear();
                        
                        // d component of trajectory
                        double df, df_dot, df_dot_dot;
                        double sf, sf_dot, sf_dot_dot;                          
                        bool use_quintic = true;
                        if (lane_change_state) {
                          //std::cout << "moving into lane " << lane_change_lane << std::endl;
                          sf = si + (((si_dot + lane_change_speed) / 2.0) * T);
                          sf_dot = lane_change_speed;
                          sf_dot_dot = 0.0;
                          df = ((double) lane_change_lane * 4.0) + 2.0;
                          df_dot = 0.0;
                          df_dot_dot = 0.0;

                          if (abs(df - di) < 0.1) {
                            std::cout << "completed lane change" << std::endl;
                            car_lane = lane_change_lane;
                            lane_change_state = false;
                          }
                        } else {
                          df = ((double) car_lane * 4.0) + 2.0;
                          df_dot = 0.0;
                          df_dot_dot = 0.0;
                          if (have_target_vehicle) {
                            double target_speed = vehicles[target_vehicle].get_velocity();
                            target_speed = (target_speed > maximum_speed) ? maximum_speed : target_speed;
                            double s_lv = vehicles[target_vehicle].get_s();             // lead vehicle initial s
                            double v_lv = vehicles[target_vehicle].get_velocity();      // lead vehicle velocity
                            double a_lv = vehicles[target_vehicle].get_acceleration();  // lead vehicle acceleration
                            double lv_fwd = (v_lv * T) + (0.5 * a_lv * pow(T,2));
                            double sf_lv = (lv_fwd > 0.0) ? (s_lv + lv_fwd) : s_lv;     // no reverse!
                            double sf_max = sf_lv - 15.0;  // don't get closer than 15m (could be made speed dependent)
                            if (sf_max <= si) {
                              std::cout << "TOO CLOSE: si=" << si << " si_dot=" << si_dot << " sf_max=" << sf_max << " veh=" <<
                                target_vehicle << " s_lv=" << s_lv << " v_lv=" << v_lv << " a_lv=" << a_lv <<
                                " sf_lv=" << sf_lv;
                              while (si >= sf_max) {
                                sf_max += 1.0;
                              }
                              std::cout << " increased sf_max to " << sf_max << std::endl;
                            }
                            assert(sf_max > si);
                            sf = sf_max;
                            sf_dot = target_speed;
                            sf_dot_dot = 0.0;
                          } else {
                            // quadric 2 eqn with no target sf
                            use_quintic = false;
                            double delta_s_dot = 0.0;  // trajectory search space
                            double target_speed = maximum_speed + delta_s_dot;
                            target_speed = (target_speed > maximum_speed) ? maximum_speed : target_speed;
                            sf_dot = target_speed;
                            sf_dot += delta_s_dot;
                            sf_dot_dot = 0.0;
                          }
                        }
                        // quintic solution for d
                        Eigen::MatrixXf m(3,3);
                        m << pow(T,3), pow(T,4), pow(T,5),
                          3.0 * pow(T,2), 4.0 * pow(T,3), 5.0 * pow(T,4),
                          6.0 * T, 12.0 * pow(T,2), 20.0 * pow(T,3);
                        Eigen::MatrixXf m_inv = m.inverse();
                        Eigen::Vector3f v;
                        v << (df - (di + (di_dot * T) + (0.5 * di_dot_dot * pow(T,2)))),
                          (df_dot - (di_dot + (di_dot_dot * T))),
                          (df_dot_dot - di_dot_dot);
                        Eigen::Vector3f soln = m_inv * v;
                        d_poly.push_back(di);
                        d_poly.push_back(di_dot);
                        d_poly.push_back(0.5 * di_dot_dot);
                        d_poly.push_back(soln(0));
                        d_poly.push_back(soln(1));
                        d_poly.push_back(soln(2));
                      
                        if (use_quintic) {
                          Eigen::Vector3f v;
                          v << (sf - (si + (si_dot * T) + (0.5 * si_dot_dot * pow(T,2)))),
                            (sf_dot - (si_dot + (si_dot_dot * T))),
                            (sf_dot_dot - si_dot_dot);
                          Eigen::Vector3f soln = m_inv * v;
                          s_poly.push_back(si);
                          s_poly.push_back(si_dot);
                          s_poly.push_back(0.5 * si_dot_dot);
                          s_poly.push_back(soln(0));
                          s_poly.push_back(soln(1));
                          s_poly.push_back(soln(2));
                        } else {
                          // quadric 2 eqn with no target sf
                          Eigen::MatrixXf m2(2,2);
                          m2 << 3.0 * pow(T,2), 4.0 * pow(T,3),
                            6.0 * T, 12.0 * pow(T,2);
                          Eigen::MatrixXf m2_inv = m2.inverse();
                          Eigen::Vector2f v2;
                          v2 << sf_dot - (si_dot + si_dot_dot * T),
                            sf_dot_dot - si_dot_dot;
                          Eigen::Vector2f soln = m2_inv * v2;
                          s_poly.push_back(si);
                          s_poly.push_back(si_dot);
                          s_poly.push_back(0.5 * si_dot_dot);
                          s_poly.push_back(soln(0));
                          s_poly.push_back(soln(1));
                          s_poly.push_back(0.0);  // use the quintic solver with alpha_5 = 0.0                        
                        }

                        //bool ok = check_path(s_poly, d_poly, sx, sy, sdx, sdy);

                        // start with copied x,y,s,d from previous solution
                        x_vals_raw = next_x_vals;
                        y_vals_raw = next_y_vals;
                        s_vals_raw = next_s_vals;
                        d_vals_raw = next_d_vals;
                        int path_length = 75;
                        assert(path_length > copy_path_cnt);
                        // keep the path size to 75 (5 copied + 70 new)
                        double s0, x0, y0;
                        max_ratio = 1.0;
                        for (int step=1; step <= path_length - copy_path_cnt; step++) {
                          auto s_soln = solve_quintic(s_poly, step * 0.02);
                          double s = get_map_s(s_soln[0]);
                          auto d_soln = solve_quintic(d_poly, step * 0.02);
                          double d = d_soln[0];
                          auto dbl_vec = getXY(s, d, sx, sy, sdx, sdy);
                          double x = dbl_vec[0];
                          double y = dbl_vec[1];

                          if (step > 1) {
                            double ratio = distance(x0, y0, x, y) / (s - s0);
                            if (ratio > max_ratio) {
                              max_ratio = ratio;
                            }
                          }
                          s0 = s;
                          x0 = x;
                          y0 = y;
                          
                          x_vals_raw.push_back(x);
                          y_vals_raw.push_back(y);
                          s_vals_raw.push_back(s);
                          d_vals_raw.push_back(d);
                        }
                        
                        bool ok = check_path(x_vals_raw, y_vals_raw);
                        if (ok) {
                          searching = false;
                        } else {
                          // try an easier trajectory
                          search_cnt++;
                          std::cout << "max_ratio=" << max_ratio << " si_dot=" << si_dot << " si_dot_dot=" << si_dot_dot
                                    << " sf_dot=" << sf_dot << " sf_dot_dot=" << sf_dot_dot << " reduce target speed "
                                    << search_cnt << " have_target "
                                    << have_target_vehicle << " lane_change " << lane_change_state << std::endl;
                          searching = false;
                        }

                      }                      
                      prev_x_vals = x_vals_raw;
                      prev_s_vals = s_vals_raw;
                      prev_d_vals = d_vals_raw;
                      next_x_vals = x_vals_raw;
                      next_y_vals = y_vals_raw;

                      if (0) {
                        std::cout << "car_x=" << car_x << " car_y=" << car_y << " next_x: ";
                        for (int i=0; i<50; i++)
                          std::cout << next_x_vals[i] << " ";
                        std::cout << std::endl;
                        std::cout << "next_y: ";
                        for (int i=0; i<50; i++)
                          std::cout << next_y_vals[i] << " ";
                        std::cout << std::endl;
                      }
                      
                      msgJson["next_x"] = next_x_vals;
                      msgJson["next_y"] = next_y_vals;

                      auto msg = "42[\"control\","+ msgJson.dump()+"]";
                      
                      //this_thread::sleep_for(chrono::milliseconds(1000));
                      ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

                      msg_cnt++;
                    }
                  } else {
                    // Manual driving
                    std::string msg = "42[\"manual\",{}]";
                    ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
                  }
                }
              });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
                    const std::string s = "<h1>Hello world!</h1>";
                    if (req.getUrl().valueLength == 1) {
                      res->end(s.data(), s.length());
                    } else {
                      // i guess this should be done more gracefully?
                      res->end(nullptr, 0);
                    }
                  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
      std::cout << "Connected!!!" << std::endl;
    });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
                      ws.close();
                      std::cout << "Disconnected" << std::endl;
                    });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}



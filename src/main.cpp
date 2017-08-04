#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <time.h>
#include <assert.h>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "Eigen-3.3/Eigen/LU"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

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
  double x = sx(s) + (d * sdx(s));
  double y = sy(s) + (d * sdy(s));
  return {x,y};
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
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

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
  
  // create a higher resolution map by interpolating the waypoint map
  std::vector<double> hires_s;
  std::vector<double> hires_x;
  std::vector<double> hires_y;
  std::vector<double> hires_dx;
  std::vector<double> hires_dy;
  for (int i=0; i < map_waypoints_s.size() - 1; i++) {
    double s0 = map_waypoints_s[i];
    double s1 = map_waypoints_s[i+1];
    int meters = s1 - s0;
    for (int j=0; j < meters; j++) {
      double s = s0 + (double) j;
      hires_s.push_back(s);
      hires_x.push_back(sx(s));
      hires_y.push_back(sy(s));
      hires_dx.push_back(sdx(s));
      hires_dy.push_back(sdy(s));
    }
  }

  std::vector<double> prev_x_vals;
  std::vector<double> prev_poly;
  int msg_cnt = 0;
  h.onMessage([&msg_cnt, &prev_x_vals, &prev_poly,
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

                      std::cout << "car_s=" << car_s << ", car_d=" << car_d << ", car_x=" << car_x << ", car_y=" << car_y << std::endl;
                      std::cout << "car_speed=" << car_speed << std::endl;

                      // Sensor Fusion Data, a list of all other cars on the same side of the road.
                      auto sensor_fusion = j[1]["sensor_fusion"];

                      json msgJson;

                      vector<double> next_x_vals;
                      vector<double> next_y_vals;
                      int copy_path_cnt = 5;
                      double car_v = mph_to_mps(car_speed);
                      double car_a = 0.0;
                      if (car_v > 0) {
                        // determine the number of points popped
                        //   then we can use prior information for the start of the next segment
                        int popped = prev_x_vals.size() - previous_path_x.size();
                        assert(popped > 0);
                        int prev_end_step = popped; // copied points are NOT counted in the polynomial step
                      
                        // stich in the beginning of the previous path and then project from there for continuity
                        int prev=0;
                        while (prev < copy_path_cnt && prev < previous_path_x.size() && previous_path_x[prev] != NULL) {
                          next_x_vals.push_back(previous_path_x[prev]);
                          next_y_vals.push_back(previous_path_y[prev]);
                          prev++;
                        }
                        if (prev > 0) {
                          assert (prev == copy_path_cnt);
                          // push only as many next points as we want to stich into the new path so the s,d values can be used
                          //   - also keep the final v,a from the previous polynomial
                          auto qsoln = solve_quintic(prev_poly, prev_end_step * 0.02);
                          car_s = qsoln[0];
                          car_v = qsoln[1];
                          car_a = qsoln[2];
                          std::cout << "start_s=" << car_s << " v=" << car_v << " a=" << car_a << std::endl;
                        }
                      } else {
                        copy_path_cnt = 0;
                      }

                      bool have_target_vehicle = false;
                      double T = 3.0;
                      double target_speed = mph_to_mps(50.0);
                      double si = car_s;
                      double si_dot = car_v;
                      double si_dot_dot = car_a;

                      double s_lv = 0.0;
                      double a_lv = 0.0;         // lead vehicle acceleration                      
                      double sf = car_s + (si_dot * T) + (0.5 * a_lv * pow(T,2));
                      // if speed is not violated, set sf to s_lv - velocity dependent buffer
                      double sf_dot = target_speed;
                      double sf_dot_dot = a_lv;

                      // std::cout << "si=" << si << ", si_dot=" << si_dot << ", si_dot_dot=" << si_dot_dot << std::endl;
                      // std::cout << "sf=" << sf << ", sf_dot=" << sf_dot << ", sf_dot_dot=" << sf_dot_dot << std::endl;
                        
                      assert(car_speed < 75.0);

                      // Using the quintic polynomial solver from Trajectory Generation
                      //   Or a quadric alternative to get the car moving suggested in the Werling paper

                      std::vector<double> poly;
                      if (have_target_vehicle) {
                      
                        Eigen::MatrixXf m(3,3);
                        m << pow(T,3), pow(T,4), pow(T,5),
                          3.0 * pow(T,2), 4.0 * pow(T,3), 5.0 * pow(T,4),
                          6.0 * T, 12.0 * pow(T,2), 20.0 * pow(T,3);
                        Eigen::MatrixXf m_inv = m.inverse();
                        Eigen::Vector3f v;
                        v << (sf - (si + si_dot * T + 0.5 * si_dot_dot * pow(T,2))),
                          (sf_dot - (si_dot + si_dot_dot * T)),
                          (sf_dot_dot - si_dot_dot);
                        Eigen::Vector3f soln = m_inv * v;
                        poly.push_back(si);
                        poly.push_back(si_dot);
                        poly.push_back(0.5 * si_dot_dot);
                        poly.push_back(soln(0));
                        poly.push_back(soln(1));
                        poly.push_back(soln(2));
                      
                      } else {
                        // quadric 2 eqn with no target sf
                        
                        double delta_s_dot = 0.0;  // trajectory search space
                        sf_dot += delta_s_dot;
                        
                        Eigen::MatrixXf m(2,2);
                        m << 3.0 * pow(T,2), 4.0 * pow(T,3),
                          6.0 * T, 12.0 * pow(T,2);
                        Eigen::MatrixXf m_inv = m.inverse();
                        Eigen::Vector2f v;
                        v << sf_dot - (si_dot + si_dot_dot * T),
                          sf_dot_dot - si_dot_dot;
                        Eigen::Vector2f soln = m_inv * v;
                        poly.push_back(si);
                        poly.push_back(si_dot);
                        poly.push_back(0.5 * si_dot_dot);
                        poly.push_back(soln(0));
                        poly.push_back(soln(1));
                        poly.push_back(0.0);  // use the quintic solver with alpha_5 = 0.0                        
                      }
                      prev_poly = poly;
                      
                      double x, y;
                      // keep the path size to 100 (10 copied + 90 new)
                      for (int step=0; step < 100 - copy_path_cnt; step++) {
                        auto qsoln = solve_quintic(poly, step * 0.02);
                        double s = qsoln[0];
                        double d = car_d;  // FIXME
                        auto dbl_vec = getXY(s, d, sx, sy, sdx, sdy);
                        x = dbl_vec[0];
                        y = dbl_vec[1];
                        next_x_vals.push_back(x);
                        next_y_vals.push_back(y);
                        //cout << " s[" << step << "] = " << s;
                      }
                      prev_x_vals = next_x_vals;
                      // std::cout << std::endl;
                      // std::cout << "next_xy: ";
                      // for (int i=0; i<next_x_vals.size(); i++) {
                      //    std::cout << "(" << next_x_vals[i] << "," << next_y_vals[i] << ") ";
                      // }
                      // std::cout << std::endl;

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



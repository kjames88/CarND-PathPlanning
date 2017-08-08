#Path Planning Project

##Path Generation

Path generation is handled in the Frenet domain.  The *s* coordinate is the longitudinal distance down the lane, and the *d* coordinate is the lateral position on the road.  Under normal circumstances, *d* should have a fixed value corresponding to the center of one of three 4m wide lanes.  The *s* coordinate has a range of up to 6945.554m, after which it wraps back to zero at the start of the track.

###Smoothing the Path

Generating x,y coordinates directly from s,d coordinates using only the provided waypoints causes high acceleration or jerk around the waypoints.  Splines are used to smooth the translation, but this a bit tricky.  Smoothing the path output from getXY does not work - the spline passes through all the points given, and the points given by getXY form a piecewise linear path.  Fitting a spline to the x,y points looks good but correlating back to s,d becomes a problem.  Some comments on the *slack* channel saved the day here:  fit functions from s to {x,y,dx,dy}.  The getXY function used here is based directly on these splines, with *dx* and *dy* required for incorporating *d*.

On the topic of *d* there may be an issue in the computation, which does not agree completely with the simulator.  While the graphics do not show a wheel on the line, the simulator periodically flags a lane violation.  *After some testing, it appears that driving 0.5m left of center performs better than driving on center.  With this bias I can clearly see the wheels on the line at times (as rendered) but there are no lane violation errors.  Conversely, when driving on center I see the wheels safely between the lines but I get lane violations.*

###Velocity without a leader

Since the car starts at zero velocity, it is immediately necessary to generate a path that gets the car moving.  Without another car in front, there is no need to specify a distance, only velocity and acceleration targets.  As suggested for this purpose in *Optimal Trajectory Generation for Dynamic Street Scenarios in a Frenet Frame* a quadric polynomial is solved for the desired maximum speed and zero acceleration at the end of a 3 second period.

###Following a leader

Using the sensor data, additional vehicles can be detected, and a close vehicle in, or entering, our lane becomes a leader vehicle.  In this case, the velocity is determined by the leader.  Target position for our car is set at a buffer distance behind the projected location of the leader at the end of the 3 second planning window.  Acceleration is not used in this location computation because negative acceleration in the simulations was often brief and projecting it out 3 seconds resulted in grossly inaccurate projections.  These projections caused strong braking and even reversal of the self driving car.  Since the path is re-computed after some tens of milliseconds, it works better to use velocity and re-compute iteratively.

When there is a leader, there is both a velocity and a position target, and for this a quintic polynomial is solved.  The quintic polynomial should minimize the jerk (derivative of acceleration).  The target *s* position is computed as the maximum of the leader's position (minus buffer distance) and our car's projected position at the average of current and target velocities.  It is very important that the target distance not be too far ahead of the vehicle, as this can cause a temporary surge in speed beyond the speed limit.

###Changing lanes

When the self driving car finds itself behind another slower-moving vehicle, it attempts to select and transition to a faster lane.  There are safety issues to consider when changing lanes; the car considers a lane blocked if it detects the following:

1. A vehicle less than 5m behind or 10m ahead (the difference due to the point location on the vehicle)
2. A vehicle projected to start from behind and end ahead of our position during the planning window
3. A vehicle starting behind and projected to be less than 5m from our projected location at the end of the planning window
4. A vehicle starting ahead and projected to be less than 20m ahead at the end of the planning window 
5. A vehicle within 50m ahead and moving slower

An important modification to observing vehicle lane positions is that a vehicle may be changing lanes.  Initially, using only the current lane (based on simple division by lane width) there were some collisions resulting from failure to account for a car that was partially over the line into our lane.  An offset is now used beyond which a vehicle is considered in both its lane and the lane it is moving toward.

Lane changes are accomplished by with a second polynomial solution, this time for *d*.  The target *d* is set to the target lane center and the *s* position is set to where our car will be at the average of current and target velocity.  Target velocity is either determined by the leading vehicle in the target lane or by a slightly reduced maximum speed if the lane is empty ahead.  Because of the lateral motion in *d*, the *s* component alone does not fully determine speed.  A precise mix is not computed; instead the nominal maximum of 50MPH is reduced by 10%.

###Trajectory Filtering

Because the trajectories are planned conservatively, I have not needed to generate multiple options and filter for jerk or acceleration.  However, because the speed limit has caused some difficulty, the paths are checked (primarily for v) and the target position is moved back if the check fails.  Attempting to reduce the target velocity was ineffective, presumably because the position target forced a velocity surge.  3 seconds of planning window limits excess acceleration issues; the cause noted for tripping this limit was braking of the leader vehicle, which when projected over 3 seconds can cause our car to brake hard or go into reverse.  (As noted above, the negative acceleration term was removed for this reason).

###Future Work

The lane change selection does not currently observe a faster or open lane two lanes away.  Thus our car will not, for example, slow down, below the speed of a vehicle in the right lane, in order to move behind a *slower* vehicle in the center lane, so that it can then move into the left lane and resume full speed travel.

# Model Predictive Control

The task of this project was to control a vehicle with MPC, instead of the PID controller used for the previous project. 

## The model
The simple kinematic model described in the lectures is used to predict the vehicle's path of motion. This model does not take forces like friction or slipping into account. The following equations are used for the state vector:
![](Images/equations.png)

## Timestep Length and Elapsed Duration
The timestep length (N) has to have enough data points to predict the path to follow accordingly. Too many points on the other side will consume too much performance. From choosing lengths between 5 and 20, I ended up using 10 points as the best option. For the elapsed duration (DT), a shorter delta time will result again in predicting a path too short. A longer delta time will produce a more "rough" path, which will possibly lead to driving in unsafe or permitted areas. I tried values between 10 and 200 ms, I determined 100 ms as the optimal value. Plus, it will come in handy for dealing with the artificial delay, that both values are the same.

## Polynomial Fitting and MPC Preprocessing
The given waypoints are transformed into the vehicle's coordinate system, to make further calculations more easy. Then these transformed points are fitted to a third degree polynomial. The coefficients can now be used to calculate the cte and epsi. They are as well used by the solver to generate the target trajectory.

## Model Predictive Control with Latency
An artificial delay of 100 ms is added to the controller code, meaning that the steering and throttle values will arrive too late, leading to inaccurate following of the optimal path as well as oscillating. As DT has been chosen to be also 100 ms, the control values can simply be shifted by one time step within the evaluator. This way, the delay for sending the control values to the vehicle can be compensated.

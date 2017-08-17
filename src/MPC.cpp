#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration
size_t N = 12;
double dt = 0.075;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

const int x0_idx = 0; // starting index of px in the vars
const int y0_idx = 1*N;
const int psi0_idx = 2*N;
const int v0_idx = 3*N;
const int cte0_idx = 4*N;
const int epsi0_idx = 5*N;
const int delta0_idx = 6*N;
const int a0_idx = 7*N - 1;

const double target_v = 65; // target speed
const double WEIGHT_V = 1.0;
const double WEIGHT_CTE = 1.0;
const double WEIGHT_PSI = 10.0;
const double WEIGHT_STEER = 2000.0;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // `vars` is a vector of variable values (state & actuators)
    // `fg` a vector of the cost fg[0] and constraints
    fg[0] = 0;

    // cost based on cte, epsi, and speed
    for(int i = 0; i < N; i++) {
      fg[0] += WEIGHT_CTE * CppAD::pow(vars[cte0_idx + i], 2);
      fg[0] += WEIGHT_PSI * CppAD::pow(vars[epsi0_idx + i], 2);
      fg[0] += WEIGHT_V * CppAD::pow(vars[v0_idx + i] - target_v, 2);
    }

    // cost based on value of actuations.
    for (int i = 0; i < N - 1; i++) {
      fg[0] += WEIGHT_STEER * CppAD::pow(vars[delta0_idx + i], 2);
      fg[0] += CppAD::pow(vars[a0_idx + i], 2);
    }

    // cost based on the value gap between sequential actuations.
    for (int i = 1; i < N - 2; i++) {
      fg[0] += CppAD::pow(vars[delta0_idx + i] - vars[delta0_idx + i -1], 2);
      fg[0] += CppAD::pow(vars[a0_idx + i] - vars[a0_idx + i -1], 2);
    }

    // Setup constrants, with initial constraints at fg[1]
    fg[1 + x0_idx] = vars[x0_idx];
    fg[1 + y0_idx] = vars[y0_idx];
    fg[1 + psi0_idx] = vars[psi0_idx];
    fg[1 + v0_idx] = vars[v0_idx];
    fg[1 + cte0_idx] = vars[cte0_idx];
    fg[1 + epsi0_idx] = vars[epsi0_idx];

    // The rest of the constraints
    for (int i = 0; i < N - 1; i++) {
      // state at time t.
      AD<double> x0 = vars[x0_idx + i];
      AD<double> y0 = vars[y0_idx + i];
      AD<double> psi0 = vars[psi0_idx + i];
      AD<double> v0 = vars[v0_idx + i];
      AD<double> cte0 = vars[cte0_idx + i];
      AD<double> epsi0 = vars[epsi0_idx + i];

      // state at time t + 1 .
      AD<double> x1 = vars[x0_idx + i + 1];
      AD<double> y1 = vars[y0_idx + i + 1];
      AD<double> psi1 = vars[psi0_idx + i + 1];
      AD<double> v1 = vars[v0_idx + i + 1];
      AD<double> cte1 = vars[cte0_idx + i + 1];
      AD<double> epsi1 = vars[epsi0_idx + i + 1];

      AD<double> delta0 = vars[delta0_idx + i];
      AD<double> a0 = vars[a0_idx + i];
      AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2]*x0*x0 + coeffs[3]*x0*x0*x0;
      AD<double> psides0 = CppAD::atan(coeffs[1]+2*coeffs[2]*x0 + 3 * coeffs[3]*x0*x0);

      fg[2 + x0_idx + i] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[2 + y0_idx + i] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[2 + psi0_idx + i] = psi1 - (psi0 + v0 * delta0 / Lf * dt);
      fg[2 + v0_idx + i] = v1 - (v0 + a0 * dt);
      fg[2 + cte0_idx + i] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
      fg[2 + epsi0_idx + i] = epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  size_t n_vars = 6 * N + 2 * (N - 1);
  // Set the number of constraints
  size_t n_constraints = 6 * N;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // Set lower and upper limits for variables. 
  // delta should be (-25DEG, +25DEG), i.e. (-0.436332, 0.436332)
  for(int i = delta0_idx; i < a0_idx; i++) {
    vars_lowerbound[i] = - 0.436332;
    vars_upperbound[i] = 0.436332;
  }

  // acceleration should (-1, 1)
  for(int i = a0_idx; i < n_vars; i++) {
    vars_lowerbound[i] = - 1.0;
    vars_upperbound[i] = 1.0;
  }

  // The rest
  for(int i = 0; i < delta0_idx; i++) {
    vars_lowerbound[i] = - 1000.0;
    vars_upperbound[i] = 1000.0;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state, which should set to the state values
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  constraints_lowerbound[x0_idx] = state[0];
  constraints_lowerbound[y0_idx] = state[1];
  constraints_lowerbound[psi0_idx] = state[2];
  constraints_lowerbound[v0_idx] = state[3];
  constraints_lowerbound[cte0_idx] = state[4];
  constraints_lowerbound[epsi0_idx] = state[5];

  constraints_upperbound[x0_idx] = state[0];
  constraints_upperbound[y0_idx] = state[1];
  constraints_upperbound[psi0_idx] = state[2];
  constraints_upperbound[v0_idx] = state[3];
  constraints_upperbound[cte0_idx] = state[4];
  constraints_upperbound[epsi0_idx] = state[5];

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // set mpc_x, mpc_y
  this->mpc_x = {};
  this->mpc_y = {};
  for (int i = 0; i < N; i++) {
    this->mpc_x.push_back(solution.x[x0_idx + i]);
    this->mpc_y.push_back(solution.x[y0_idx + i]);
  }

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  vector<double> result = {solution.x[delta0_idx], solution.x[a0_idx]};

  return result;
}

Eigen::VectorXd MPC::GetState(Eigen::VectorXd coeffs, double v, double steering_angle, double throttle)
{
  // Vehicle model
  Eigen::VectorXd result = Eigen::VectorXd(6);

  // cte is the deviation from center y
  // after local coord transformation, vehicle is at (0, 0)
  // to account for latency, predicate state ahead by 0.1 second
  double latency = 0.1;

  double px = 0.0, py = 0.0;
  double cte = (coeffs[0]+coeffs[1]*px+coeffs[2]*px*px+coeffs[3]*px*px*px) - py;
  double epsi = - atan(coeffs[1]);
  double px_pred = v * latency;
  double py_pred = 0;
  double psi_pred = - v * steering_angle * latency / Lf;
  double v_pred = v + throttle * latency;
  double cte_pred = cte + v * sin(epsi) * latency;
  double epsi_pred = epsi + psi_pred;

  result << px_pred, py_pred, psi_pred, v_pred, cte_pred, epsi_pred;
  // cout << "State: " << result << endl;
  return result;
}

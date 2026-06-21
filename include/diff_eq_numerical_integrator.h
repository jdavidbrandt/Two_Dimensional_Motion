#ifndef DIFF_EQ_NUMERICAL_INTEGRATOR_H_INCLUDED
#define DIFF_EQ_NUMERICAL_INTEGRATOR_H_INCLUDED

#include <cmath>
#include <vector>
#include <iostream>
#include <functional>
#include <cstdarg>
#include <stdexcept>
#include <algorithm>
#include <tuple>


//InitialCondition(double independent, double dependent)
//take two doubles, one being the independent, the second being the dependent, and assigns them to pair. y(0) = 1 would construct this tuple ->     (0, 1)
class InitialCondition{
    public:
        
    // pairs is the tuple (independent, dependent)
        std::tuple<double, double> pairs;
        double independent;
        double dependent;

        //constructor, assigns above members with what is passed in via initializer
        InitialCondition(double passed_independent_0, double passed_dependent_0)
            :pairs(passed_independent_0, passed_dependent_0),
            independent(passed_independent_0),
            dependent(passed_dependent_0)    
        {


        }

};

//base class
class DiffEqn {
    public:
        //each differential equation will have a data set, the predicted independent, and dependent vals
        std::vector<double> predicted_rk4_independent_vals;
        std::vector<double> predicted_rk4_dependent_vals;
};

//FirstOrderODE(std::function<double(double,double)> passed_equation, InitialCondition initial_condition)
//of form dy/dt + P(t)y = Q(t). Define Diff eqn by solving for dy/dt algebraically, and returning what is on the right hand side of the eqn in main. dy/dt = return value
//this is done in main()
class FirstOrderODE: public DiffEqn{
    public:
        //define initial conditions, initiate a function (empty, it is defined in main), and initiate vectors for independent and dependent predicted values
        std::function<double(double,double)> equation;
        InitialCondition initialCondition;

        //FirstOrderODE(std::function<double(double,double)> passed_equation, InitialCondition initial_condition)
        //takes in what dy/dt is equal to (passed_equation), and an object of type InitialCondition, and assigns passed_equation, and initial conditions via initializer
        FirstOrderODE(std::function<double(double,double)> passed_equation, InitialCondition initial_condition)
            : equation(passed_equation),
            initialCondition(initial_condition.independent, initial_condition.dependent)

        {
            //sets the initial condition
            predicted_rk4_independent_vals.push_back(initial_condition.independent);
            predicted_rk4_dependent_vals.push_back(initial_condition.dependent);         
        }

    //Add first order Euler Method
        double FirstOrderEulerMethodGetPos(double step_size, InitialCondition initial_condition){
            double forward_pos = initial_condition.dependent + step_size*equation(initial_condition.independent, initial_condition.dependent);
            return forward_pos;
        }
    //FirstOrderRK4Solve(double step_size, double independent_final)
    //function to solve the the First Order Linear Diff Eqn, given a step size, and independent final starting from t = 0
        void FirstOrderRK4Solve(double step_size, double independent_final){

            //Runge-Kutta constants initialization, and setting the first step to the initial conditions
            double k1, k2, k3, k4;
            double independent_step = initialCondition.independent;
            double dependent_step = initialCondition.dependent;
            

            //Define number of steps, and start estimating slopes at every step, give them a weighted average, and fill vectors so that one can plot predicted values
            int n = (independent_final - initialCondition.independent)/step_size;
            for (int i=0; i<n; ++i){

                //if we ever divide by 0, go to the next step, and don't calculate k vals
                if(std::isnan(equation(independent_step, dependent_step)) || std::isinf(equation(independent_step, dependent_step))){
                    independent_step+=step_size;
                }
                else{
                    k1 = step_size * equation(independent_step, dependent_step);
                    k2 = step_size * equation(independent_step + 0.5*step_size, dependent_step + 0.5*k1);
                    k3 = step_size * equation(independent_step + 0.5*step_size, dependent_step + 0.5*k2);
                    k4 = step_size * equation(independent_step + step_size, dependent_step + k3);
                    independent_step += step_size;
                    dependent_step += (k1 + 2*k2 + 2*k3 + k4)/6.0;

                    //adds each predicted point to two vectors, such that one can plot.
                    predicted_rk4_independent_vals.push_back(independent_step);
                    predicted_rk4_dependent_vals.push_back(dependent_step);
                }   
            }
        }
        //function overload for choosing an independent plotting start value that's not the "initial condition".
        //The initial condition must lie on the curve, it doesn't have to be the start of the curve
        void FirstOrderRK4Solve(double step_size, double independent_plot_start_val, double independent_final){

            //Runge-Kutta constants initialization, and setting the first step to the initial conditions
            double l1, l2, l3, l4;
            double independent_step = initialCondition.independent;
            double dependent_step = initialCondition.dependent;

            //Define number of steps, and start estimating slopes at every step, give them a weighted average, and fill vectors so that one can plot predicted values
            int m = (initialCondition.independent - independent_plot_start_val)/step_size;
            for(int j = 0; j < m; ++j){
                //if we ever divide by 0, go to the next step, and don't calculate l vals
                if(std::isnan(equation(independent_step, dependent_step)) || std::isinf(equation(independent_step, dependent_step))){
                    independent_step-=step_size;
                }
                else{
                    l1 = step_size * equation(independent_step, dependent_step);
                    l2 = step_size * equation(independent_step + 0.5*step_size, dependent_step + 0.5*l1);
                    l3 = step_size * equation(independent_step + 0.5*step_size, dependent_step + 0.5*l2);
                    l4 = step_size * equation(independent_step + step_size, dependent_step + l3);

                    //increment
                    independent_step -= step_size;
                    dependent_step -= (l1 + 2*l2 + 2*l3 + l4)/6.0;

                    //store, such that one can plot.
                    predicted_rk4_independent_vals.push_back(independent_step);
                    predicted_rk4_dependent_vals.push_back(dependent_step);
                }
            }
            
            //since I am decrementing, i must reverse the order of predicted_rk4_*_vals so that the first and middle points are not connected, such that I can plot left to right, instead of the middle out.
            std::reverse(predicted_rk4_independent_vals.begin(), predicted_rk4_independent_vals.end());
            std::reverse(predicted_rk4_dependent_vals.begin(), predicted_rk4_dependent_vals.end());

            //call non-overloaded function to plot the back portion of the set
            FirstOrderRK4Solve(step_size, independent_final);
        }
};

//SecondOrderODE(std::function<double(double,double,double)> passed_equation, InitialCondition posCondition, InitialCondition velCondition)
//of form y'' = f(t, y, y')
class SecondOrderODE: public DiffEqn{
    public:
        //define initial conditions, initiate equations used, and initiate two additional vectors for y' values, predicted_prime_vals
        std::function<double(double,double,double)> v_dot; 
        InitialCondition posInitialCondition;
        InitialCondition velInitialCondition;
        std::vector<double> predicted_prime_rk4_independent_vals, predicted_prime_rk4_dependent_vals;
        
        //Takes in an equation, y''=f(t,y,y'), and two InitialCondition, one for position, the other for velocity. the initializes all via initializer
        SecondOrderODE(std::function<double(double,double,double)> passed_equation, InitialCondition posCondition, InitialCondition velCondition)
            : v_dot(passed_equation),
            posInitialCondition(posCondition.independent, posCondition.dependent),
            velInitialCondition(velCondition.independent, velCondition.dependent)
            
        {
            
            //The equations are broken into two first order, coupled differential equations upon invoking SecondOrderRK4Solve() where v = f1(t, y, y_dot) and v_dot = f2(t, y, y_dot) 
            //sets the first steps
            predicted_rk4_independent_vals.push_back(posCondition.independent);
            predicted_rk4_dependent_vals.push_back(posCondition.dependent);
            predicted_prime_rk4_independent_vals.push_back(velCondition.independent);
            predicted_prime_rk4_dependent_vals.push_back(velCondition.dependent);


            
        }
        
        //SecondOrderRK4Solve(double step_size, double independent_final)
        //takes in step_size, and independent_final, and populates members with doubles such that one can plot
        void SecondOrderRK4Solve(double step_size, double independent_final){
            //initialize Runge-Kutta constants k being for dependent, and l being for dependent_prime, along with the substitution equation, v
            std::function<double(double,double,double)> v = [](double t, double y, double y_dot){
                return y_dot;
            };
            double k1, k2, k3, k4, l1, l2, l3, l4;
            double independent_step = posInitialCondition.independent;
            double dependent_step = posInitialCondition.dependent;
            double independent_prime_step = velInitialCondition.independent;
            double dependent_prime_step = velInitialCondition.dependent;


            //Define number of steps, and start estimating slopes at every step, give them a weighted average, and fill vectors so that one can plot predicted values
            int n = (independent_final - posInitialCondition.independent)/step_size;
            for (int i = 0; i < n; ++i){
                //if we ever divide by 0, go to the next step, and don't calculate k and l vals
                if(std::isnan(v_dot(independent_step, dependent_step, dependent_prime_step)) || std::isinf(v_dot(independent_step, dependent_step, dependent_prime_step))){
                    independent_step+=step_size;
                    independent_prime_step+=step_size;
                }
                else{
                    k1 = step_size*v(independent_step, dependent_step, dependent_prime_step);
                    l1 = step_size*v_dot(independent_step, dependent_step, dependent_prime_step);
                    k2 = step_size*v(independent_step + 0.5*step_size, dependent_step + 0.5*step_size*k1, dependent_prime_step + 0.5*step_size*l1);
                    l2 = step_size*v_dot(independent_step + 0.5*step_size, dependent_step + 0.5*step_size*k1, dependent_prime_step + 0.5*step_size*l1);
                    k3 = step_size*v(independent_step + 0.5*step_size, dependent_step + 0.5*step_size*k2, dependent_prime_step + 0.5*step_size*l2);
                    l3 = step_size*v_dot(independent_step + 0.5*step_size, dependent_step + 0.5*step_size*k2, dependent_prime_step + 0.5*step_size*l2);
                    k4 = step_size*v(independent_step + step_size, dependent_step + step_size*k3, dependent_prime_step + step_size*l3);
                    l4 = step_size*v_dot(independent_step + step_size, dependent_step + step_size*k3, dependent_prime_step + step_size*l3);
                    
                

                    //increment
                    //independent_step += 0.001;
                    //independent_prime_step += step_size;

                    independent_step = posInitialCondition.independent + (i + 1)*step_size;
                    independent_prime_step = velInitialCondition.independent + (i + 1)*step_size;

                    dependent_step += (k1 + 2*k2 + 2*k3 + k4)/6.0;
                    dependent_prime_step += (l1 + 2*l2 + 2*l3 + l4)/6.0;

                    //store, such that one can plot
                    predicted_rk4_independent_vals.push_back(independent_step);
                    predicted_rk4_dependent_vals.push_back(dependent_step);
                    predicted_prime_rk4_independent_vals.push_back(independent_prime_step);
                    predicted_prime_rk4_dependent_vals.push_back(dependent_prime_step);
                }
            }
        }

        //takes step size, independent_initial (from initial conditions), and independent_plot_starting_val, (where the plot actually starts plotting)
        //function overload to pass one more arg, then call the non-overloaded function to complete the plot from initial conditions to independent_final
        void SecondOrderRK4Solve(double step_size, double independent_plot_starting_val, double independent_final){
            //substitution variable which is always true, let v = y_dot    
            std::function<double(double,double,double)> v = [](double t, double y, double y_dot){
                return y_dot;
            };
            double i1, i2, i3, i4, j1, j2, j3, j4;

            //have them start at the initial conditions, and work backwards to the independent_plot_starting_val
            double independent_step = posInitialCondition.independent;
            double dependent_step = posInitialCondition.dependent;
            double independent_prime_step = velInitialCondition.independent;
            double dependent_prime_step = velInitialCondition.dependent;

            //Define number of steps, and start estimating slopes at every step, give them a weighted average, and fill vectors so that one can plot predicted values
            int m = (posInitialCondition.independent - independent_plot_starting_val)/step_size;
            for(int j = 0; j < m; ++j){
                //if we ever divide by 0, go to the next step, and don't calculate i and j vals
                if(std::isnan(v_dot(independent_step, dependent_step, dependent_prime_step)) || std::isinf(v_dot(independent_step, dependent_step, dependent_prime_step))){
                    independent_step-=step_size;
                    independent_prime_step-=step_size;
                }
                else{
                    i1 = step_size*v(independent_step, dependent_step, dependent_prime_step);
                    j1 = step_size*v_dot(independent_step, dependent_step, dependent_prime_step);
                    i2 = step_size*v(independent_step + 0.5*step_size, dependent_step + 0.5*step_size*i1, dependent_prime_step + 0.5*step_size*j1);
                    j2 = step_size*v_dot(independent_step + 0.5*step_size, dependent_step + 0.5*step_size*i1, dependent_prime_step + 0.5*step_size*j1);
                    i3 = step_size*v(independent_step + 0.5*step_size, dependent_step + 0.5*step_size*i2, dependent_prime_step + 0.5*step_size*j2);
                    j3 = step_size*v_dot(independent_step + 0.5*step_size, dependent_step + 0.5*step_size*i2, dependent_prime_step + 0.5*step_size*j2);
                    i4 = step_size*v(independent_step + step_size, dependent_step + step_size*i3, dependent_prime_step + step_size*j3);
                    j4 = step_size*v_dot(independent_step + step_size, dependent_step + step_size*i3, dependent_prime_step + step_size*j3);

                    //decrement
                    independent_step -= step_size;
                    independent_prime_step -= step_size;
                    dependent_step -= (i1 + 2*i2 + 2*i3 + i4)/6.0;
                    dependent_prime_step -= (j1 + 2*j2 + 2*j3 + j4)/6.0;

                    //storing such that one can plot
                    predicted_rk4_independent_vals.push_back(independent_step);
                    predicted_rk4_dependent_vals.push_back(dependent_step);
                    predicted_prime_rk4_independent_vals.push_back(independent_prime_step);
                    predicted_prime_rk4_dependent_vals.push_back(dependent_prime_step);
                }
            }
            
            //since I am decrementing, i must reverse the order of predicted_rk4_*_vals so that the first and middle points are not connect, such that I can plot left to right, instead of the middle out.
            std::reverse(predicted_rk4_independent_vals.begin(), predicted_rk4_independent_vals.end());
            std::reverse(predicted_rk4_dependent_vals.begin(), predicted_rk4_dependent_vals.end());
            std::reverse(predicted_prime_rk4_independent_vals.begin(), predicted_prime_rk4_independent_vals.end());
            std::reverse(predicted_prime_rk4_dependent_vals.begin(), predicted_prime_rk4_dependent_vals.end());

            //call non-overloaded function to plot the back half of the set
            SecondOrderRK4Solve(step_size, independent_final);
        }
        
            
};

//Create a new class called SecondOrderCoupledODE that has two equations. Probably don't do inheritance as it would be confusing which equation is for which thing
class CoupledSecondOrderODE{
    public:
        //x'' = f(t, x, y, x', y') and y'' = g(t, x, y, x', y')
        //By convention, let's use f for x'' equation, and g for y'' equation
        std::function<double(double,double,double,double,double)> vx_dot;
        std::function<double(double,double,double,double,double)> vy_dot;
        InitialCondition f_PositionIC;
        InitialCondition f_VelocityIC;
        InitialCondition g_PositionIC;
        InitialCondition g_VelocityIC;

        std::vector<double> f_independent_vals;
        std::vector<double> f_dependent_vals;
        std::vector<double> f_prime_independent_vals;
        std::vector<double> f_prime_dependent_vals;
        std::vector<double> g_independent_vals;
        std::vector<double> g_dependent_vals;
        std::vector<double> g_prime_independent_vals;
        std::vector<double> g_prime_dependent_vals;

        CoupledSecondOrderODE(std::function<double(double,double,double,double,double)> f,
                              std::function<double(double,double,double,double,double)> g,
                              InitialCondition f_PositionIC,
                              InitialCondition f_VelocityIC,
                              InitialCondition g_PositionIC,
                              InitialCondition g_VelocityIC)
                            : vx_dot(f),
                            vy_dot(g),
                            f_PositionIC(f_PositionIC),
                            f_VelocityIC(f_VelocityIC),
                            g_PositionIC(g_PositionIC),
                            g_VelocityIC(g_VelocityIC)
        {
            f_independent_vals.push_back(f_PositionIC.independent);
            f_dependent_vals.push_back(f_PositionIC.dependent);
            f_prime_independent_vals.push_back(f_VelocityIC.independent);
            f_prime_dependent_vals.push_back(f_VelocityIC.dependent);
            g_independent_vals.push_back(g_PositionIC.independent);
            g_dependent_vals.push_back(g_PositionIC.dependent);
            g_prime_independent_vals.push_back(g_VelocityIC.independent);
            g_prime_dependent_vals.push_back(g_VelocityIC.dependent);

        };

        void CoupledSecondOrderRK4Solve(double step_size, double independent_final){
            std::function<double(double,double,double,double,double)> vx = [](double t, double x, double y, double x_dot, double y_dot){
                return x_dot;
            };

            std::function<double(double,double,double,double,double)> vy = [](double t, double x, double y, double x_dot, double y_dot){
                return y_dot;
            };

            double k1, k2, k3, k4, l1, l2, l3, l4, m1, m2, m3, m4, n1, n2, n3, n4;
            double f_independent_step = f_PositionIC.independent;
            double f_dependent_step = f_PositionIC.dependent;
            double f_prime_independent_step = f_VelocityIC.independent;
            double f_prime_dependent_step = f_VelocityIC.dependent;

            double g_independent_step = g_PositionIC.independent;
            double g_dependent_step = g_PositionIC.dependent;
            double g_prime_independent_step = g_VelocityIC.independent;
            double g_prime_dependent_step = g_VelocityIC.dependent;

            int N = (independent_final - f_PositionIC.independent) / step_size;
            
            for(int i = 0; i < N; ++i){
                //if we ever divide by 0, go to the next step and don't plot
                if(std::isnan(vx_dot(f_independent_step, f_dependent_step, g_dependent_step,f_prime_dependent_step, g_prime_dependent_step)) || std::isinf(vx_dot(f_independent_step, f_dependent_step, g_dependent_step,f_prime_dependent_step, g_prime_dependent_step))){
                    f_independent_step += step_size;
                    f_prime_independent_step += step_size;
                    g_independent_step += step_size;
                    g_prime_independent_step += step_size;
                }

                if(std::isnan(vy_dot(f_independent_step, f_dependent_step, g_dependent_step,f_prime_dependent_step, g_prime_dependent_step)) || std::isinf(vy_dot(f_independent_step, f_dependent_step, g_dependent_step,f_prime_dependent_step, g_prime_dependent_step))){
                    f_independent_step += step_size;
                    f_prime_independent_step += step_size;
                    g_independent_step += step_size;
                    g_prime_independent_step += step_size;
                }
                


                k1 = vx(f_independent_step, f_dependent_step, g_dependent_step, f_prime_dependent_step, g_prime_dependent_step)*step_size;
                k2 = vy(g_independent_step, f_dependent_step, g_dependent_step, f_prime_dependent_step, g_prime_dependent_step)*step_size;
                k3 = vx_dot(f_independent_step, f_dependent_step, g_dependent_step, f_prime_dependent_step, g_prime_dependent_step)*step_size;
                k4 = vy_dot(f_independent_step, f_dependent_step, g_dependent_step, f_prime_dependent_step, g_prime_dependent_step)*step_size;

                l1 = (vx(f_independent_step, f_dependent_step, g_dependent_step, f_prime_dependent_step, g_prime_dependent_step) + k3)*step_size;
                l2 = (vy(f_independent_step, f_dependent_step, g_dependent_step, f_prime_dependent_step, g_prime_dependent_step) + k4)*step_size;
                l3 = vx_dot(f_independent_step + 0.5*step_size, f_dependent_step + 0.5*k1, g_dependent_step + 0.5*k2, f_prime_dependent_step + 0.5*k3, g_prime_dependent_step + 0.5*k4)*step_size;
                l4 = vy_dot(f_independent_step + 0.5*step_size, f_dependent_step + 0.5*k1, g_dependent_step + 0.5*k2, f_prime_dependent_step + 0.5*k3, g_prime_dependent_step + 0.5*k4)*step_size;

                m1 = (vx(f_independent_step, f_dependent_step, g_dependent_step, f_prime_dependent_step, g_prime_dependent_step) + l3)*step_size;
                m2 = (vy(f_independent_step, f_dependent_step, g_dependent_step, f_prime_dependent_step, g_prime_dependent_step) + l4)*step_size;
                m3 = vx_dot(f_independent_step + 0.5*step_size, f_dependent_step + 0.5*l1, g_dependent_step + 0.5*l2, f_prime_dependent_step + 0.5*l3, g_prime_dependent_step + 0.5*l4)*step_size;
                m4 = vy_dot(f_independent_step + 0.5*step_size, f_dependent_step + 0.5*l1, g_dependent_step + 0.5*l2, f_prime_dependent_step + 0.5*l3, g_prime_dependent_step + 0.5*l4)*step_size;

                n1 = (vx(f_independent_step, f_dependent_step, g_dependent_step, f_prime_dependent_step, g_prime_dependent_step) + m3)*step_size;
                n2 = (vy(f_independent_step, f_dependent_step, g_dependent_step, f_prime_dependent_step, g_prime_dependent_step) + m4)*step_size;
                n3 = vx_dot(f_independent_step + 0.5*step_size, f_dependent_step + 0.5*m1, g_dependent_step + 0.5*m2, f_prime_dependent_step + 0.5*m3, g_prime_dependent_step + 0.5*m4)*step_size;
                n4 = vy_dot(f_independent_step + 0.5*step_size, f_dependent_step + 0.5*m1, g_dependent_step + 0.5*m2, f_prime_dependent_step + 0.5*m3, g_prime_dependent_step + 0.5*m4)*step_size;

                //increment
                f_independent_step += step_size;
                f_prime_independent_step += step_size;
                g_independent_step += step_size;
                g_prime_independent_step += step_size;

                //Weighted Average
                f_dependent_step += (k1 + 2*l1 + 2* m1 + n1)/6.0;
                g_dependent_step += (k2 + 2*l2 + 2*m2 + n2)/6.0;
                f_prime_dependent_step += (k3 + 2*l3 + 2*m3 + n3)/6.0;
                g_prime_dependent_step += (k4 + 2*l4 + 2*m4 + n4)/6.0;

                f_independent_vals.push_back(f_independent_step);
                f_dependent_vals.push_back(f_dependent_step);
                f_prime_independent_vals.push_back(f_prime_independent_step);
                f_prime_dependent_vals.push_back(f_prime_dependent_step);

                g_independent_vals.push_back(g_independent_step);
                g_dependent_vals.push_back(g_dependent_step);
                g_prime_independent_vals.push_back(g_prime_independent_step);
                g_prime_dependent_vals.push_back(g_prime_dependent_step);



            }
        }


};


////////
//In Progress: Not Functional
//Non-functional
class HigherOrderODE: public DiffEqn{
    public:
        //takes the order of the differential equation, the passed equation (with coefficent of 1 for highest order), and a vector of initial conditions sorted from 0th derivative to nth order derivative
        int order;
        std::vector<InitialCondition> vInitialConditions;

        HigherOrderODE(int passed_order, std::function<double(double,double)> passed_equation, std::vector<InitialCondition> passed_InitialConditions){
            order = passed_order;
            //equation = passed_equation;
            vInitialConditions = passed_InitialConditions;

        }
    
};








#endif
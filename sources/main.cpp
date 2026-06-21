#include <iostream>
#include "raylib.h"
#include "../include/diff_eq_numerical_integrator.h"
#include "../include/physics_objects.h"
#include "rlgl.h"
#include "raymath.h"
#include <cmath>
#include <vector>
#include <matplot/matplot.h>


#define SCREEN_WIDTH (1920)
#define SCREEN_HEIGHT (1080)
#define WINDOW_TITLE "Two Dimensional Motion"
#define FPS 60





//SampleAtTime(double t, const std::vector<double>& times, const std::vector<double>& values)
//Takes in a time to be sampled at t, times, and values, and grabs the time t and returns the value
double SampleAtTime(double t, const std::vector<double>& times, const std::vector<double>& values){
    //safeguard for if times or values happens to be empty
    if(times.size() < 2) return values.empty() ? 0.0 : values.front();
    if(t <= times.front()) return values.front();
    if(t >= times.back()) return values.back();


    const double h = times[1] - times[0];                           //equivalent to step size
    const double fidx = (t - times.front()) / h;                    //the amount of samples between time t and the initial time 
    const size_t i = static_cast<size_t>(fidx);                     
    const double frac = fidx - static_cast<double>(i);

    return values[i] + frac * (values[i+1] - values[i]);

    
}

//Important to note, 1 meter == 1 pixel
int main()
{
    const double pi = 3.1415926535;
    InitWindow(SCREEN_WIDTH, SCREEN_HEIGHT, WINDOW_TITLE);
    SetTargetFPS(FPS);

    //Circle properties
    double blue_diameter = 10e-5;   //meters
    double red_diameter = 10.0;      //meters
    double circ_xpos = (SCREEN_WIDTH / 2.0) - (blue_diameter / 2.0);
    double circ_ypos = (SCREEN_HEIGHT / 2.0) + (blue_diameter);
    double circ_mass = 1.0;       //kgs
    double red_circ_mass = 1000.0; //kgs

    //Create circle for this one since there's air, don't use rectangle
    //We are going to use the same position for both blue and red circles, but not draw the blue one due to limitations
    PhysicsCircle blue_circle = {circ_xpos, circ_ypos, blue_diameter / 2.0, BLUE, circ_mass};
    PhysicsCircle red_circle = {circ_xpos, circ_ypos, red_diameter / 2.0, RED, red_circ_mass};


    const double G = 6.67430e-11;                                                   //Newton's Gravitational Constant
    const double mass_of_planet = 5.9722e24;                                       //Earth
    const double radius_of_planet = 6378137;                                      //meters (equatorial radius for simplicity)
    double gravity = -1*(G * mass_of_planet) / (pow(radius_of_planet, 2));

    //Camera properties, set to center of the screen
    Vector2 offset = {(SCREEN_WIDTH / 2.0), (SCREEN_HEIGHT / 2.0)};
    Vector2 target = {(SCREEN_WIDTH / 2.0), (SCREEN_HEIGHT / 2.0)};
    float rotation = 0.0f;

    //zooming in, because 1 pixel = 1 meter
    float zoom = 1.0f;
    Camera2D camera = {offset, target, rotation, zoom};

    
    
    //Take into account the ratio of the quadratic drag force compared to the linear drag force (at STP) (Classical Mechanics, Taylor, pg 45)
    //f_quad / f_lin = γDv/ β   ->   (1.6 * 10^3)*Dv, <- (AT STP) where D is the diameter of the sphere, and v is the velocity at which it is moving
    //For a baseball of D = 7 cm, with v = 5 m/s,                       f_quad / f_lin ~ 600,   so f_lin can be ignored,
    //For a drop of oil, of D = 1.5 μm, and v = 5e-5 m/s,               f_quad / f_lin ~ 1e_7,  so f_quad can be ignored
    //For a raindrop of D = 1 mm and v = 0.6 m/s,                       f_quad / f_lin ~ 1,     so neither can be neglected

    //BLUE CIRCLE: FLIN used, FQUAD ignored
    //For the blue ball of diameter 10e-5 m and velocity 10e-5 m/s
    //Unfortunately becuase of the way raylib draws by the pixel, zooming in to much smaller than a pixel causes the animations to be clunky
    //But we can still plot this information to compare to the quadratic air resistance, and we will


    //RED CIRCLE: FQUAD used, FLIN ignored
    //For the blue ball, of diameter 10m and velocity starting at 0 and getting up to even 5 m/s
    //f_quad / f_lin ~ 40,000 even going at only 5 m/s, we can ignore f_lin here

    //ORANGE CIRCLE: Both FQUAD and FLIN used

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Define differential equation of motion here, (the y-axis in raylib is flipped, we will take care of this later)    
    //Linear Air Resistance Only, is a first order ODE, because the forces do not depend on the position, only the velocity v
    // x'' = f(t, x, x')
    std::function<double(double,double,double)> xBlue_Equation = [blue_circle](double t, double x, double x_dot){
        return (-1*blue_circle.b*x_dot) / blue_circle.mass;
    };

    // y'' = g(t, y, y')
    std::function<double(double,double,double)> yBlue_Equation = [blue_circle, gravity](double t, double y, double y_dot){
        return(gravity - (blue_circle.b * y_dot) / blue_circle.mass);
    };

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Quadratic Air Resistance Only
    // x'' = f(t, x, y, x', y')
    std::function<double(double,double,double,double, double)> xRed_Equation = [red_circle, gravity](double t, double x, double y, double x_dot, double y_dot){
        return -1*red_circle.c*(sqrt(x_dot*x_dot + y_dot*y_dot)*x_dot) / red_circle.mass;
    };

    std::function<double(double,double,double,double,double)> yRed_Equation = [red_circle, gravity](double t, double x, double y, double x_dot, double y_dot){
        return gravity - ((red_circle.c*(sqrt(x_dot*x_dot + y_dot*y_dot)))*y_dot)/red_circle.mass;
    };


    //Initial conditions and values, generic t_initial and t_final for all rectangles
    double t_initial = 0.0;
    double t_final = 25;
    double step_size = 0.001;


    //Blue Circle
    //y(0) = blue_circly.y,  y'(0) = sin(pi/4) and x(0) = blue_circle.x, x'(0) = cos(pi/4)
    //Create Initial condition objects here
    double blue_xinitial = blue_circle.x;
    double blue_vxinitial = 10e-5*cos(pi/4.0);
    double blue_yinitial = blue_circle.y;
    double blue_vyinitial = 10e-5*sin(pi/4);
    InitialCondition xPosCondition(t_initial, blue_xinitial);
    InitialCondition xVelCondition(t_initial, blue_vxinitial);
    InitialCondition yPosCondition(t_initial, blue_yinitial);
    InitialCondition yVelCondition(t_initial, blue_vyinitial);
   
    SecondOrderODE xBlue_ODE(xBlue_Equation, xPosCondition, xVelCondition);
    SecondOrderODE yBlue_ODE(yBlue_Equation, yPosCondition, yVelCondition);

    xBlue_ODE.SecondOrderRK4Solve(step_size, t_final);
    yBlue_ODE.SecondOrderRK4Solve(step_size, t_final);
     
    //SecondOrderODE blue_falling_body(blue_eqn_of_motion, initial_condition_gravity_pos, initial_condition_gravity_vel);
    //blue_falling_body.SecondOrderRK4Solve(step_size, t_final);
    std::vector<double> blue_xtime_vals = xBlue_ODE.predicted_rk4_independent_vals;
    std::vector<double> blue_xpos_vals = xBlue_ODE.predicted_rk4_dependent_vals;
    std::vector<double> blue_ytime_vals = yBlue_ODE.predicted_rk4_independent_vals;
    std::vector<double> blue_ypos_vals = yBlue_ODE.predicted_rk4_dependent_vals;
    //std::vector<double> vel_time_vals = blue_falling_body.predicted_prime_rk4_independent_vals;
    //std::vector<double> vel_vals = blue_falling_body.predicted_prime_rk4_dependent_vals;

    //Red Circle, Quadratic, must use CoupledSecondOrder to solve
    double red_xinitial = red_circle.x;
    double red_yinitial = red_circle.y;
    double red_vxinitial = 50;
    double red_vyinitial = 50;
    InitialCondition Red_xPosCondition(t_initial, red_xinitial);
    InitialCondition Red_xVelCondition(t_initial, red_vxinitial);
    InitialCondition Red_yPosCondition(t_initial, red_yinitial);
    InitialCondition Red_yVelCondition(t_initial, red_vyinitial);

    CoupledSecondOrderODE Red_CoupledODE(xRed_Equation, yRed_Equation, Red_xPosCondition, Red_xVelCondition, Red_yPosCondition, Red_yVelCondition);
    Red_CoupledODE.CoupledSecondOrderRK4Solve(step_size, t_final);

    std::vector<double> red_xtime_vals = Red_CoupledODE.f_independent_vals;
    std::vector<double> red_xpos_vals = Red_CoupledODE.f_dependent_vals;
    std::vector<double> red_ytime_vals = Red_CoupledODE.g_independent_vals;
    std::vector<double> red_ypos_vals = Red_CoupledODE.g_dependent_vals;
    std::vector<double> red_vxtime_vals = Red_CoupledODE.f_prime_independent_vals;
    std::vector<double> red_vytime_vals = Red_CoupledODE.g_prime_independent_vals;
    std::vector<double> red_vxvel_vals = Red_CoupledODE.f_prime_dependent_vals;
    std::vector<double> red_vyvel_vals = Red_CoupledODE.g_prime_dependent_vals;



    double elapsed = 0.0;
    //const double rect_start_y = initial_condition_gravity_pos.independent;
    
     const int x_offset = 10;

    while (!WindowShouldClose())
    {
        BeginDrawing();
        ClearBackground(BLACK);

        BeginMode2D(camera);
            
        //Can't draw due to limitations with Raylib ):

        //Diff Eqn Blue Rect, flipping y because intrinsically, +y is down, so I'm flipping it such that -y is down
        //const double blue_start_y = (SCREEN_HEIGHT / 2.0) - (rect_height / 2.0);
        //blue_rect_ypos = SampleAtTime(elapsed, time_vals, y_pos_vals);
        //double bworld_x = SampleAtTime(elapsed, xtime_vals, xpos_vals);
        //double screen_x = x_initial + (x_initial - world_x);
        //double bworld_y = SampleAtTime(elapsed, ytime_vals, ypos_vals);
        //double bscreen_y = blue_yinitial + (blue_yinitial - bworld_y);
        //DrawRectangle(rect_xpos + x_offset, screen_y, rect_width, rect_height, BLUE);
        //DrawCircle(bworld_x, bscreen_y, blue_diameter/2.0, BLUE);

        
        //Quad
        const double red_start_y = (SCREEN_HEIGHT / 2.0) - (red_diameter / 2.0);
        double red_worldx = SampleAtTime(elapsed, red_xtime_vals, red_xpos_vals);
        double red_worldy = SampleAtTime(elapsed, red_ytime_vals, red_ypos_vals);
        double red_screeny = red_yinitial +(red_yinitial - red_worldy);

        DrawCircle(red_worldx, red_screeny, red_circle.radius, RED);


        
        elapsed += GetFrameTime();



        EndMode2D();
        EndDrawing();
        if(elapsed > 26){
            break;
        }
        
    }


    CloseWindow();

//Testing Coupled Second Order Diff EQN
    
    // Plotting RK4 predictions
    auto fig = matplot::figure(true);
    matplot::plot(red_xtime_vals, red_xpos_vals, "r",
                    red_ytime_vals, red_ypos_vals, "r");
    
    
        

    std::vector<std::string> labels = {"RK4-Blue", "RK4-Orange", "Actual" };
    
    
    matplot::xlabel("Time (s)");
    matplot::ylabel("Height (m)");
    matplot::title("Free Falling Body, Δ = 0.00001");
    //matplot::legend(labels);

    fig->width(800);
    fig->height(600);

    //matplot::save("plots/delta.png");

    matplot::show(fig);
    
    return 0;
}



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



//Important to note, 1 meter := 1 pixel
int main()
{
    //Define Constants, and Camera Properties
    const double pi = 3.1415926535;
    const double G = 6.67430e-11;                                                   //Newton's Gravitational Constant
    const double mass_of_planet = 5.9722e24;                                       //Earth
    const double radius_of_planet = 6378137;                                      //meters (equatorial radius for simplicity)
    double gravity = -1*(G * mass_of_planet) / (pow(radius_of_planet, 2));       // 9.81 m/s^2
    Vector2 offset = {(SCREEN_WIDTH / 2.0), (SCREEN_HEIGHT / 2.0)};
    Vector2 target = {(SCREEN_WIDTH / 2.0), (SCREEN_HEIGHT / 2.0)};
    float rotation = 0.0f;
    float zoom = 10.0f;
    Camera2D camera = {offset, target, rotation, zoom};

    InitWindow(SCREEN_WIDTH, SCREEN_HEIGHT, WINDOW_TITLE);
    SetTargetFPS(FPS);

    //Circle properties, all circles start from center, and we will give them an x offset when drawn
    double blue_diameter = 1;       //meters
    double red_diameter = 1;       //meters
    double orange_diameter = 1;   //meters
    double circ_xpos = 0.0;
    double circ_ypos = 0.0;
    double blue_circ_mass = 1;              //kgs
    double red_circ_mass = 10;             //kgs
    double orange_circ_mass = 1000.0;     //kgs

    //We are going to use the same position for both blue and red circles, but not draw the blue one due to limitations with raylib
    PhysicsCircle blue_circle = {circ_xpos, circ_ypos, blue_diameter / 2.0, BLUE, blue_circ_mass};
    PhysicsCircle red_circle = {circ_xpos, circ_ypos, red_diameter / 2.0, RED, red_circ_mass};
    PhysicsCircle orange_circle = {circ_xpos, circ_ypos, orange_diameter / 2.0, ORANGE, orange_circ_mass};
    
    //Notice I am using a circle, this is important, because the differential equations would be different for a different shape
    //Also note that since we are dealing with non-conservative forces (drag), we cannot use Lagrangian Mechanics
    //Take into account the ratio of the quadratic drag force compared to the linear drag force (at STP) (Classical Mechanics, Taylor, pg 45)
    //f_quad / f_lin = γDv/ β   ->   (1.6 * 10^3)*Dv, <- (AT STP) where D is the diameter of the sphere, and v is the velocity at which it is moving
    //For a baseball of D = 7 cm, with v = 5 m/s,                       f_quad / f_lin ~ 600,   so f_lin can be ignored,
    //For a drop of oil, of D = 1.5 μm, and v = 5e-5 m/s,               f_quad / f_lin ~ 1e_7,  so f_quad can be ignored
    //For a raindrop of D = 1 mm and v = 0.6 m/s,                       f_quad / f_lin ~ 1,     so neither can be neglected

    //Typically, at STP, the only time linear resistance comes into play is if you start from rest, and the medium is viscous, so what is mostly studied here is quadratic drag.
    //Since we are starting at ~50 m/s, let's compare masses here with quadratic drag, and look at terminal velocity, comparing masses and diameters
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // IMPORTANT NOTE: I AM LEAVING THESE SEPARATELY DEFINED, BECAUSE I AM GOING TO COME BACK AND COMPARE LINEAR DRAG, SO I WILL BE CHANGING ONE OR TWO OF THESE EQUATIONS TO INCLUDE LINEAR AIR RESISTANCE
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Define differential equations of motion for circles here, (the y-axis in raylib is flipped, we will take care of this when drawing the circle later)  
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //BLUE CIRCLE: 
    // x'' = f(t, x, y, x', y'), and y''(t) = g(t, x, y, x', y')
    std::function<double(double,double,double,double,double)> xBlue_Equation = [blue_circle](double t, double x, double y, double x_dot, double y_dot){
        return (-1 * ((blue_circle.c)*sqrt(x_dot*x_dot + y_dot*y_dot)*x_dot) / blue_circle.mass);
    };

    std::function<double(double,double,double,double,double)> yBlue_Equation = [blue_circle, gravity](double t, double x, double y, double x_dot, double y_dot){
        return (gravity - ((blue_circle.c* sqrt(x_dot*x_dot + y_dot*y_dot) * y_dot) / blue_circle.mass));
    };

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //RED CIRCLE: 
    std::function<double(double,double,double,double, double)> xRed_Equation = [red_circle, gravity](double t, double x, double y, double x_dot, double y_dot){
        return  -1*red_circle.c*(sqrt(x_dot*x_dot + y_dot*y_dot)*x_dot) / red_circle.mass;
    };

    std::function<double(double,double,double,double,double)> yRed_Equation = [red_circle, gravity](double t, double x, double y, double x_dot, double y_dot){
        return gravity - ((red_circle.c*(sqrt(x_dot*x_dot + y_dot*y_dot)))*y_dot)/red_circle.mass;
    };
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //ORANGE CIRCLE: 
    std::function<double(double,double,double,double,double)> xOrange_Equation = [orange_circle, gravity](double t, double x, double y, double x_dot, double y_dot){
        return (-1 * ((orange_circle.c)*sqrt(x_dot*x_dot + y_dot*y_dot)*x_dot) / orange_circle.mass);
    };

    std::function<double(double,double,double,double,double)> yOrange_Equation = [orange_circle, gravity](double t, double x, double y, double x_dot, double y_dot){
        return (gravity - ((orange_circle.c* sqrt(x_dot*x_dot + y_dot*y_dot) * y_dot) / orange_circle.mass));
    };


    //Initial conditions and values, generic t_initial and t_final for all circles
    double t_initial = 0.0;
    double t_final = 25;
    double step_size = 0.001;

    //Blue Circle, initial conditions are starting at origin (0,0) with 50 x̂ + 50 ŷ which gives |v| ~= 70.71
    double blue_xinitial = blue_circle.x;
    double blue_vxinitial = 50;
    double blue_yinitial = blue_circle.y;
    double blue_vyinitial = 50;
    InitialCondition Blue_xPosCondition(t_initial, blue_xinitial);
    InitialCondition Blue_xVelCondition(t_initial, blue_vxinitial);
    InitialCondition Blue_yPosCondition(t_initial, blue_yinitial);
    InitialCondition Blue_yVelCondition(t_initial, blue_vyinitial);
   
    //Create the coupled differential equations, and solve them usin RK 4th Order
    CoupledSecondOrderODE Blue_CoupledODE(xBlue_Equation, yBlue_Equation, Blue_xPosCondition, Blue_xVelCondition, Blue_yPosCondition, Blue_yVelCondition);
    Blue_CoupledODE.CoupledSecondOrderRK4Solve(step_size, t_final);
    
    //define plot vals where f is associated with the x axis, and g is associated with the y axis
    std::vector<double> blue_xtime_vals = Blue_CoupledODE.f_independent_vals;
    std::vector<double> blue_xpos_vals = Blue_CoupledODE.f_dependent_vals;
    std::vector<double> blue_ytime_vals = Blue_CoupledODE.g_independent_vals;
    std::vector<double> blue_ypos_vals = Blue_CoupledODE.g_dependent_vals;
    std::vector<double> blue_vxtime_vals = Blue_CoupledODE.f_prime_independent_vals;
    std::vector<double> blue_vytime_vals = Blue_CoupledODE.g_prime_independent_vals;
    std::vector<double> blue_vxvel_vals = Blue_CoupledODE.f_prime_dependent_vals;
    std::vector<double> blue_vyvel_vals = Blue_CoupledODE.g_prime_dependent_vals;

    //Red Circle, repeated above
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

    //Orange Circle, repeated above
    double orange_xinitial = orange_circle.x;
    double orange_yinitial = orange_circle.y;
    double orange_vxinitial = 50;
    double orange_vyinitial = 50;
    InitialCondition Orange_xPosCondition(t_initial, orange_xinitial);
    InitialCondition Orange_xVelCondition(t_initial, orange_vxinitial);
    InitialCondition Orange_yPosCondition(t_initial, orange_yinitial);
    InitialCondition Orange_yVelCondition(t_initial, orange_vyinitial);

    CoupledSecondOrderODE Orange_CoupledODE(xOrange_Equation, yOrange_Equation, Orange_xPosCondition, Orange_xVelCondition, Orange_yPosCondition, Orange_yVelCondition);
    Orange_CoupledODE.CoupledSecondOrderRK4Solve(step_size, t_final);

    std::vector<double> orange_xtime_vals = Orange_CoupledODE.f_independent_vals;
    std::vector<double> orange_xpos_vals = Orange_CoupledODE.f_dependent_vals;
    std::vector<double> orange_ytime_vals = Orange_CoupledODE.g_independent_vals;
    std::vector<double> orange_ypos_vals = Orange_CoupledODE.g_dependent_vals;
    std::vector<double> orange_vxtime_vals = Orange_CoupledODE.f_prime_independent_vals;
    std::vector<double> orange_vytime_vals = Orange_CoupledODE.g_prime_independent_vals;
    std::vector<double> orange_vxvel_vals = Orange_CoupledODE.f_prime_dependent_vals;
    std::vector<double> orange_vyvel_vals = Orange_CoupledODE.g_prime_dependent_vals;
    


    //Track time
    double elapsed = 0.0;
    double y_offset = 540;
    double x_offset = 955;
    
    while (!WindowShouldClose())
    {
        BeginDrawing();
        ClearBackground(BLACK);

        BeginMode2D(camera);
                    
        //Define coordinates and draw circles at every sampled time
        double blue_worldx = SampleAtTime(elapsed, blue_xtime_vals, blue_xpos_vals);
        double blue_worldy = SampleAtTime(elapsed, blue_ytime_vals, blue_ypos_vals);
        double blue_screeny = blue_yinitial + (blue_yinitial - blue_worldy);

        double red_worldx = SampleAtTime(elapsed, red_xtime_vals, red_xpos_vals);
        double red_worldy = SampleAtTime(elapsed, red_ytime_vals, red_ypos_vals);
        double red_screeny = red_yinitial +(red_yinitial - red_worldy);

        double orange_worldx = SampleAtTime(elapsed, orange_xtime_vals, orange_xpos_vals);
        double orange_worldy = SampleAtTime(elapsed, orange_ytime_vals, orange_ypos_vals);
        double orange_screeny = orange_yinitial +(orange_yinitial - orange_worldy);

        //Creating vectors to pass into DrawCircleV()
        Vector2 blue_vector(blue_worldx + x_offset, blue_screeny + y_offset);
        Vector2 red_vector(red_worldx + x_offset, red_screeny + y_offset);
        Vector2 orange_vector(orange_worldx + x_offset, orange_screeny + y_offset);

        //Switched to DrawCircleV() which takes in vectors to eliminate int clunky-ness from DrawCircle()
        DrawCircleV(blue_vector, blue_circle.radius, BLUE);
        DrawCircleV(red_vector, red_circle.radius, RED);
        DrawCircleV(orange_vector, orange_circle.radius, ORANGE);

        //Add time
        elapsed += GetFrameTime();

        EndMode2D();
        EndDrawing();

        if(elapsed > t_final){
            break;
        }
        
    }


    CloseWindow();

//Testing Coupled Second Order Diff EQN
    
    // Plotting RK4 predictions

    ////////////////////////////////////////////////////////////////////////
    //POSITION//////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    auto posfig = matplot::figure(true);
    auto posaxes = matplot::axes();
    matplot::plot(blue_ytime_vals, blue_ypos_vals, "b", 
        blue_xtime_vals, blue_xpos_vals, "b--",
        red_ytime_vals, red_ypos_vals, "r",
        red_xtime_vals, red_xpos_vals, "r--",
        orange_ytime_vals, orange_ypos_vals, "k",
        orange_xtime_vals, orange_xpos_vals, "k--");
    posaxes->xlim({t_initial, t_final});
    std::vector<std::string> pos_labels = {"Blue Y", "Blue X", "Red Y", "Red X", "Orange Y", "Orange X"};
    //std::vector<std::string> pos_labels = {"Red Y", "Red X", "Orange Y", "Orange X"};
    auto pos_legend = matplot::legend(pos_labels);
    pos_legend->location(matplot::legend::general_alignment::bottomleft);
    matplot::xlabel("Time (s)");
    matplot::ylabel("Position (m)");
    matplot::title("Position, Quadratic Air Resistance, Varying Diameter");
    //matplot::save("plots/Pos_Quad_AR_VaryingDiameterNoBlue.png");
    matplot::show(posfig);

    /////////////////////////////////////////////////////////////////////////
    //VELOCITY//////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    auto velfig = matplot::figure(true);
    auto velaxes = matplot::axes();
    matplot::plot(blue_vytime_vals, blue_vyvel_vals, "b", 
        blue_vxtime_vals, blue_vxvel_vals, "b--",
        red_vytime_vals, red_vyvel_vals, "r",
        red_vxtime_vals, red_vxvel_vals, "r--",
        orange_vytime_vals, orange_vyvel_vals, "k",
        orange_vxtime_vals, orange_vxvel_vals, "k--");
    velaxes->xlim({t_initial, t_final});
    std::vector<std::string> vel_labels = {"Blue Y", "Blue X", "Red Y", "Red X", "Orange Y", "Orange X"};
    //std::vector<std::string> vel_labels = {"Red Y", "Red X", "Orange Y", "Orange X"};
    auto vel_legend = matplot::legend(vel_labels);
    vel_legend->location(matplot::legend::general_alignment::right);
    matplot::xlabel("Time (s)");
    matplot::ylabel("Velocity (m/s)");
    matplot::title("Velocity, Quadratic Air Resistance, Varying Diameter");
    //matplot::save("plots/Vel_Quad_AR_VaryingDiameterNoBlue.png");
    matplot::show(velfig);

    ////////////////////////////////////////////////////////////////
    //LEAVE PLOTTING BELOW FOR LINEAR REVISIT///////////////////////
    ////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////
    //BLUE Circle///////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////

    // auto bluefig = matplot::figure(true);
    // auto blueaxes = matplot::axes();
    // matplot::plot(blue_ytime_vals, blue_ypos_vals, "b", 
    //     blue_xtime_vals, blue_xpos_vals, "b--");
    // blueaxes->xlim({t_initial, t_final});
    // std::vector<std::string> blue_labels = {"Y Position", "X Position"};
    // auto blue_legend = matplot::legend(blue_labels);
    // blue_legend->location(matplot::legend::general_alignment::bottomleft);
    // matplot::xlabel("Time (s)");
    // matplot::ylabel("Position (m)");
    // matplot::title("Blue Position, Air Resistance, m = 0.1");
    // //matplot::save("plots/delta.png");
    // matplot::show(bluefig);

    // auto bluefig1 = matplot::figure(true);
    // auto blueaxes1 = matplot::axes();
    // matplot::plot(blue_vytime_vals, blue_vyvel_vals, "b",
    //     blue_vxtime_vals, blue_vxvel_vals, "b--");
    // blueaxes1->xlim({t_initial, t_final});
    // std::vector<std::string> blue_labels1 = {"Y Veloctiy", "X Velocity"};
    // auto blue_legend1 = matplot::legend(blue_labels1);
    // blue_legend1->location(matplot::legend::general_alignment::right);
    // matplot::xlabel("Time (s)");
    // matplot::ylabel("Velocity (m/s)");
    // matplot::title("Blue Velocity, Air Resistance, m = 0.001");
    // matplot::show(bluefig1);


    ////////////////////////////////////////////////////////////////
    //RED Circle////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////
    // auto redfig = matplot::figure(true);
    // auto redaxes = matplot::axes();
    // matplot::plot(red_ytime_vals, red_ypos_vals, "r", 
    //     red_xtime_vals, red_xpos_vals, "r--");
    // redaxes->xlim({t_initial, t_final});
    // std::vector<std::string> red_labels = {"Y Position", "X Position"};
    // matplot::legend(red_labels);
    // auto red_legend = matplot::legend(red_labels);
    // red_legend->location(matplot::legend::general_alignment::bottomleft);
    // matplot::xlabel("Time (s)");
    // matplot::ylabel("Position (m)");
    // matplot::title("Red Position, Air Resistance, m = 1");
    // //matplot::save("plots/delta.png");
    // matplot::show(redfig);

    // auto redfig1 = matplot::figure(true);
    // auto redaxes1 = matplot::axes();
    // matplot::plot(red_vytime_vals, red_vyvel_vals, "r",
    //     red_vxtime_vals, red_vxvel_vals, "r--");
    // redaxes1->xlim({t_initial, t_final});
    // std::vector<std::string> red_labels1 = {"Y Veloctiy", "X Velocity"};
    // auto red_legend1 = matplot::legend(red_labels1);
    // matplot::xlabel("Time (s)");
    // matplot::ylabel("Velocity (m/s)");
    // matplot::title("Red Velocity, Air Resistance, m = 1");
    // matplot::show(redfig1);


    //////////////////////////////////////////////////////////////////
    //ORANGE Circle///////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////
    // auto orangefig = matplot::figure(true);
    // auto orangeaxes = matplot::axes();
    // matplot::plot(orange_ytime_vals, orange_ypos_vals, "k", 
    //     orange_xtime_vals, orange_xpos_vals, "k--");
    // orangeaxes->xlim({t_initial, t_final});
    // std::vector<std::string> orange_labels = {"Y Position", "X Position"};
    // auto orange_legend = matplot::legend(orange_labels);
    // orange_legend->location(matplot::legend::general_alignment::right);
    // matplot::xlabel("Time (s)");
    // matplot::ylabel("Position (m)");
    // matplot::title("Orange Position, Air Resistance, m = 1000");
    // //matplot::save("plots/delta.png");
    // matplot::show(orangefig);

    // auto orangefig1 = matplot::figure(true);
    // auto orangeaxes1 = matplot::axes();
    // matplot::plot(orange_vytime_vals, orange_vyvel_vals, "k",
    //     orange_vxtime_vals, orange_vxvel_vals, "k--");
    // orangeaxes1->xlim({t_initial, t_final});
    // std::vector<std::string> orange_labels1 = {"Y Veloctiy", "X Velocity"};
    // auto orange_legend1 = matplot::legend(orange_labels1);
    // orange_legend1->location(matplot::legend::general_alignment::right);
    // matplot::xlabel("Time (s)");
    // matplot::ylabel("Velocity (m/s)");
    // matplot::title("Orange Velocity, Air Resistance, m = 1000");
    // matplot::show(orangefig1);

    return 0;
}



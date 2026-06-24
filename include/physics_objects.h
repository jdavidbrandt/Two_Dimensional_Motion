#ifndef PHYSICS_OBJECTS_H_INCLUDED
#define PHYSICS_OBJECTS_H_INCLUDED
#include <raylib.h>

inline constexpr double kBetaLin = 1.6e-4;                  // linear:    f_lin  = b*v,  b = β·D
inline constexpr double kGammaQuad = 0.25;                  // quadratic: f_quad = c*v², c = γ·D²


struct PhysicsRectangle : public Rectangle{
    float mass;

    //override of rectangle that assigns mass
    PhysicsRectangle(float x, float y, float width, float height, float mass){
        this->x = x;
        this->y = y;
        this->width = width;
        this->height = height;
        this->mass = mass;
    }
};

struct PhysicsCircle{
    double mass;
    double x;
    double y;
    double radius;
    Color color;
    double b;
    double c;

    PhysicsCircle(double centerX, double centerY, double radius, Color color, double mass)
        : x(centerX),
        y(centerY),
        radius(radius),
        color(color),
        mass(mass)

    {

        b = kBetaLin * 2*radius;                            // linear:    f_lin  = b*v,  b = β·D
        c = kGammaQuad * pow(2*radius,2);                   // quadratic: f_quad = c*v², c = γ·D²
    }
    
};



#endif
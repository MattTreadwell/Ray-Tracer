//
// hw3
//
//     Created: 4/20/20
//      Author: Matthew Treadwell
// 
//

#ifndef HW3_HW3_H
#define HW3_HW3_H

#ifdef WIN32
#include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
#include <GL/gl.h>
#include <GL/glut.h>
#elif defined(__APPLE__)
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> // for sqrt
#include <algorithm> // for max,min
#include <limits>
#include <stdexcept>
#include <omp.h>

#ifdef WIN32
#define strcasecmp _stricmp
#endif

#include <imageIO.h>

/*
 * Macros
 */

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0
#define DEG_TO_RAD (M_PI / 180.0)

// Minimum distance that two objects must be from each other (including itself) aka epsilon
#define effective_zero  1e-10
// Square bounds for SSAA (0.5 == 1 pixel width/height)
#define AA_BOUND 0.5
// Optimal grid rotation angle for SSAA
#define RGSS_ANGLE 26.6

/*
 * Provided Structs
 */

// Represent objects triangles, spheres, and lights from scene description (provided)
struct Vertex {
    double position[3];
    double color_diffuse[3];
    double color_specular[3];
    double normal[3];
    double shininess;
};

struct Triangle {
    Vertex v[3];
};

struct Sphere {
    double position[3];
    double color_diffuse[3];
    double color_specular[3];
    double shininess;
    double radius;
};

struct Light {
    double position[3];
    double color[3];
};

/*
 * User-defined structs
 */
enum objType {
    SPHERE, TRIANGLE, ERROR
};

// For a given ray, we will generate a clipped color (almost identical to vector)
struct Color {
    double r, g, b;

    Color() {}

    Color(double _r, double _g, double _b) {
        // Clamp to range 0.0-1.0
        r = clamp(_r);
        g = clamp(_g);
        b = clamp(_b);
    }

    // addition operator to add color components together
    Color operator+(const Color &rhs) const {
        // Ensure values are always clamped by returning a new color object
        return Color(r + rhs.r, g + rhs.g, b + rhs.b);
    }

    Color operator-(const Color &rhs) const {
        return Color(r - rhs.r, g - rhs.g, b - rhs.b);
    }

    Color operator*(double scalar) const {
        return Color(r * scalar, g * scalar, b * scalar);
    }

    Color operator/(double div) const {
        return Color(r / div, g / div, b / div);
    }

    static double clamp(double val) {
        // clamp to range 0.0-1.0
        return std::max(0.0, std::min(val, 1.0));
    }
};

// https://joncraton.org/blog/67/simple-vector-library-c/ (for class definition)
struct Vector {
    double x, y, z;

    Vector() {}

    Vector(double _x, double _y, double _z)
            : x(_x), y(_y), z(_z) {}

    Vector(const Vertex& v)
            : x(v.position[0]), y(v.position[1]), z(v.position[2]) {}

    Vector(const Vertex& v, bool normalVector)
            : x(v.normal[0]), y(v.normal[1]), z(v.normal[2]) {}

    Vector(const Light& l)
            : x(l.position[0]), y(l.position[1]), z(l.position[2]) {}

    Vector(const Sphere& s)
            : x(s.position[0]), y(s.position[1]), z(s.position[2]) {}

    // These operators will greatly clean up the math involved
    // Using const references will speed up performance
    bool operator==(const Vector &rhs) const {
        return (rhs.x == x) && (rhs.y == y) && (rhs.z == z);
    }

    Vector operator+(const Vector &rhs) const {
        return Vector(x + rhs.x, y + rhs.y, z + rhs.z);
    }

    Vector operator-(const Vector &rhs) const {
        return Vector(x - rhs.x, y - rhs.y, z - rhs.z);
    }

    Vector operator*(double scalar) const {
        return Vector(x * scalar, y * scalar, z * scalar);
    }

    Vector operator/(double scalar) const {
        return Vector(x / scalar, y / scalar, z / scalar);
    }

    Vector cross(const Vector &rhs) const {
        return Vector(y * rhs.z - z * rhs.y,
                      z * rhs.x - x * rhs.z,
                      x * rhs.y - y * rhs.x);
    }

    double dot(const Vector &rhs) const {
        return (x * rhs.x) + (y * rhs.y) + (z * rhs.z);
    }

    double length() const {
        return sqrt(x * x + y * y + z * z);
    }

    void normalize() {
        // Normalize the vector
        double magnitude = length();

        x /= magnitude;
        y /= magnitude;
        z /= magnitude;
    }
};

/*
 * user-defined class definitions
 */

class Phong {
public:
    // static methods that give the phong shading values for a given light, object
    static Color genShade(const Sphere &sphere, const Light &light, const Vector &intersect);

    static Color genShade(const Triangle &triangle, const Light &light, const Vector &intersect);
};

class Ray {
private:
    // Initial defined values
    Vector direction, origin;

    // vars generated from genIntersect
    Vector min_intersect;
    size_t min_index; // array index of the closest intersecting object
    objType intersectType = ERROR;

    // Get distance of closest intersection, returns negative if intersection does not exist
    double getIntersectionDistance(const Sphere &sphere) const;

    double getIntersectionDistance(const Triangle &triangle) const;

public:

    Ray() {}

    Ray(const Vector &_origin, const Vector &_direction)
            : origin(_origin), direction(_direction) {}

    // Intersection Code
    /*** These functions have no side-effects ***/
    bool intersectExists() const {
        return intersectType != ERROR;
    }

    // Depends on genIntersect to be called first
    Color fireShadowRays(const Sphere *spheres, const Triangle *triangles, const Light *lights, size_t num_spheres,
                         size_t num_triangles, size_t num_lights) const;

    /*** These functions have side-effects (stored in members variables) ***/
    // Modifies: intersectType, min_distance, min_index
    void genIntersect(const Sphere *spheres, const Triangle *triangles, size_t num_spheres, size_t num_triangles);
};


/*
 * Globals
 * Note: All globals used in raytracing must be const in order for algorithm to be parallelized
 */

// Note that the tan function in math.h takes RADIANS as argument, not degrees.
const double aspect_ratio = double(WIDTH) / double(HEIGHT);
const double tan_fov = tan((fov / 2) * DEG_TO_RAD);

char *filename = NULL;

// Represents the final color for each pixel that will be output
unsigned char buffer[HEIGHT][WIDTH][3];

// Arrays that hold all of the scene objects
Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

/*
 * Function Prototypes
 */

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b);

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b);

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b);


#endif //HW3_HW3_H

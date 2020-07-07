/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Matthew Treadwell
 * *************************
*/

#include "raytrace.h"
/*
 * User-defined class functions
 */

Color Phong::genShade(const Sphere &sphere, const Light &light, const Vector &intersect) {
    Color retVal;

    Vector center(sphere);
    Vector lightOrigin(light);

    Vector L = lightOrigin - intersect;
    L.normalize();
    Vector N = intersect - center;
    N.normalize();
    double L_dot_N = Color::clamp(L.dot(N));

    Vector R = N*(2*L_dot_N) - L;
    R.normalize();
    Vector V = intersect * -1.0;
    V.normalize();
    double R_dot_V = Color::clamp(R.dot(V));

    // I = lightColor * (kd * (L dot N) + ks * (R dot V) ^ sh)
    // kd == obj.color_diffuse[x]
    // ks == obj.color_specular[x]
    // sh == obj.shininess
    // if L dot N < 0, you should clamp L dot N to zero; same for R dot V
    double r = light.color[0] * (sphere.color_diffuse[0] * L_dot_N + sphere.color_specular[0] * pow(R_dot_V, sphere.shininess));
    double g = light.color[1] * (sphere.color_diffuse[1] * L_dot_N + sphere.color_specular[1] * pow(R_dot_V, sphere.shininess));
    double b = light.color[2] * (sphere.color_diffuse[2] * L_dot_N + sphere.color_specular[2] * pow(R_dot_V, sphere.shininess));

    retVal = Color(r,g,b);

    return retVal;
}

Color Phong::genShade(const Triangle &triangle, const Light &light, const Vector &intersect) {
    // Reference: Lecture 16, starting at slide 16
    Color retVal;

    Vector A = Vector(triangle.v[0]);
    Vector B = Vector(triangle.v[1]);
    Vector C = Vector(triangle.v[2]);
    Vector n = (B-A).cross(C-A);

    // Barycentric coords
    double alpha = n.dot((C - B).cross(intersect - B)) / n.dot(n);
    double beta = n.dot((A - C).cross(intersect - C)) / n.dot(n);
    double gamma = 1.0 - alpha - beta;

    Vector lightOrigin(light);
    Vector L = lightOrigin - intersect;
    L.normalize();
    Vector N = Vector(triangle.v[0],true)*alpha + Vector(triangle.v[1],true)*beta + Vector(triangle.v[2],true)*gamma;
    N.normalize();
    double L_dot_N = Color::clamp(L.dot(N));

    Vector R = N*(2*L_dot_N) - L;
    R.normalize();
    Vector V = intersect * -1.0;
    V.normalize();
    double R_dot_V = Color::clamp(R.dot(V));

    // Need to interpolate triangle components
    double color_diffuse[3];
    double color_specular[3];

    for(size_t i = 0; i < 3; ++i) {
        color_diffuse[i] = triangle.v[0].color_diffuse[i]*alpha
                + triangle.v[1].color_diffuse[i]*beta
                + triangle.v[2].color_diffuse[i]*gamma;
        color_specular[i] = triangle.v[0].color_specular[i]*alpha
                + triangle.v[1].color_specular[i]*beta
                + triangle.v[2].color_specular[i]*gamma;
    }

    double shininess = triangle.v[0].shininess*alpha + triangle.v[1].shininess*beta + triangle.v[2].shininess*gamma;

    // I = lightColor * (kd * (L dot N) + ks * (R dot V) ^ sh)
    // kd == color_diffuse[x]
    // ks == color_specular[x]
    // sh == shininess
    double r = light.color[0] * (color_diffuse[0] * L_dot_N + color_specular[0] * pow(R_dot_V, shininess));
    double g = light.color[1] * (color_diffuse[1] * L_dot_N + color_specular[1] * pow(R_dot_V, shininess));
    double b = light.color[2] * (color_diffuse[2] * L_dot_N + color_specular[2] * pow(R_dot_V, shininess));

    retVal = Color(r,g,b);

    return retVal;
}

double Ray::getIntersectionDistance(const Sphere &sphere) const {
    // Reference: Lecture 16, starting at slide 5

    Vector center(sphere.position[0], sphere.position[1], sphere.position[2]);
    Vector dist = origin - center;
    double r = sphere.radius;
    double b = 2 * direction.dot(dist);
    double c = dist.dot(dist) - pow(r, 2);
    double disc = pow(b, 2) - 4 * c;

    // Imaginary number
    if(disc < 0) {
        return -1.0;
    }

    // Calculate both quadratic results
    double t0 = (-b + sqrt(disc)) / 2.0;
    double t1 = (-b - sqrt(disc)) / 2.0;
    double t;

    if (t0 < effective_zero && t1 < effective_zero) {
        return -1.0;
    } else if (t0 < effective_zero) {
        t = t1;
    } else if (t1 < effective_zero) {
        t = t0;
    } else {
        t = std::min(t0, t1);
    }

    return t;
}

double Ray::getIntersectionDistance(const Triangle &triangle) const {
    // Reference: Lecture 16, starting at slide 10
    Vector A = Vector(triangle.v[0]);
    Vector B = Vector(triangle.v[1]);
    Vector C = Vector(triangle.v[2]);

    Vector n = (B - A).cross(C - A);
    n.normalize();
    double d = -n.dot(A);
    double t = -(n.dot(origin) + d) / n.dot(direction);
    Vector intersect = origin + direction * t;

    // Create 3 triangles and check if they are parallel
    double t1 = n.dot((C - B).cross(intersect - B));
    double t2 = n.dot((A - C).cross(intersect - C));
    double t3 = n.dot((B - A).cross(intersect - A));

    // Intersection is invalid if:
    // 1. If n dot d = 0, no intersection (ray parallel to plane)
    // 2. If t leq 0, the intersection is behind ray origin
    // 3. fails point-in-triangle test for any of the three triangles (slide 19 & 20)
    if (abs(n.dot(direction)) <= effective_zero ||      // 1
        t < effective_zero ||                           // 2
        t1 < 0.0 || t2 < 0.0 || t3 < 0.0) {
        return -1.0;
    }

    return t;
}

Color Ray::fireShadowRays(const Sphere *spheres, const Triangle *triangles, const Light *lights, size_t num_spheres,
                          size_t num_triangles, size_t num_lights) const {
    Color retVal(0.0,0.0,0.0);

    if(!intersectExists()) {
        throw std::invalid_argument("Error: fireShadowRays called without an intersection");
    }

    // For each light
    Vector lightOrigin, lightDirection;
    for(size_t i = 0; i < num_lights; ++i) {
        lightOrigin = Vector(lights[i].position[0], lights[i].position[1], lights[i].position[2]);
        lightDirection = Vector(lightOrigin - min_intersect);
        lightDirection.normalize();
        Ray shadowRay(min_intersect, lightDirection);

        bool blocked = false;

        // Get closest intersecting object with shadow ray
        shadowRay.genIntersect(spheres, triangles, num_spheres, num_triangles);

        // If there is an intersection, check if the intersection is before the light
        if(shadowRay.intersectExists()) {
            double shadow_len = (shadowRay.min_intersect - min_intersect).length();
            double light_len = (lightOrigin - min_intersect).length();
            // if the intersecting object is before the light, the light is blocked
            if(shadow_len < light_len) {
                blocked = true;
            }
        }

        // if this shadow ray is unblocked, evaluate local Phong model
        if(!blocked) {
            switch (intersectType) {
                case SPHERE:
                    retVal = retVal + Phong::genShade(spheres[min_index], lights[i], min_intersect);
                    break;
                case TRIANGLE:
                    retVal = retVal + Phong::genShade(triangles[min_index], lights[i], min_intersect);
                    break;
                default:
                    throw std::invalid_argument("Error: fireShadowRays called without an intersection");
            }
        }
    }

    return retVal;
}

void Ray::genIntersect(const Sphere *spheres, const Triangle *triangles, size_t num_spheres, size_t num_triangles) {
    // Finds the closest intersecting object for this given ray, if it exists
    min_index = 0;
    intersectType = ERROR;
    double min_distance = std::numeric_limits<double>::max();
    double tmp_distance = std::numeric_limits<double>::max();

    // Spheres
    for (size_t i = 0; i < num_spheres; ++i) {
        // If intersection
        tmp_distance = getIntersectionDistance(spheres[i]);
        if (tmp_distance > 0.0 && tmp_distance < min_distance) {
            intersectType = SPHERE;
            min_distance = tmp_distance;
            min_index = i;
            min_intersect = origin + direction * min_distance;
        }
    }

    // Triangles
    for (size_t i = 0; i < num_triangles; ++i) {
        // If intersection
        tmp_distance = getIntersectionDistance(triangles[i]);
        if (tmp_distance > 0.0 && tmp_distance < min_distance) {
            intersectType = TRIANGLE;
            min_distance = tmp_distance;
            min_index = i;
            min_intersect = origin + direction * min_distance;
        }
    }
}

/*
 * Color generation code
 */

// Generates color for a given ray
Color genColor(Ray ray) {
    Color retVal(0,0,0);

    // Generate phong component
    ray.genIntersect(spheres, triangles, num_spheres, num_triangles);
    if(ray.intersectExists()) {
        Color phong = ray.fireShadowRays(spheres, triangles, lights, num_spheres, num_triangles, num_lights);
        Color reflection (0.0,0.0,0.0);

        // Add the phong and reflection components
        retVal = retVal + phong;
        retVal = retVal + reflection;
    } else {
        // background color
        retVal = Color(1.0,1.0,1.0);
    }

    return retVal;
}

// Generate a ray for a given pixel
Ray genRay(double x, double y) {
    Ray retVal;

    double screen_x = 2 * ((x + 0.5) / WIDTH) - 1;
    double screen_y = 2 * ((y + 0.5) / HEIGHT) - 1;

    Vector direction = Vector(screen_x * aspect_ratio * tan_fov, screen_y * tan_fov, -1);
    direction.normalize();
    Vector origin = Vector(0.0, 0.0, 0.0);

    retVal = Ray(origin, direction);

    return retVal;
}

/*
 * Start of user code
 */
void draw_scene() {
    // Render every column in parallel
    #pragma omp parallel for
    for (unsigned int x = 0; x < WIDTH; x++) {
        double deg;
        for (unsigned int y = 0; y < HEIGHT; y++) {
            /**** For each pixel x,y ****/
            // https://en.wikipedia.org/wiki/Supersampling
            // 4x SSAA using rotated grid
            // For an optimal pattern, the rotation angle is arctan(1/2) or 26.6 degrees
            Color anti_aliased(0.0,0.0,0.0);
            for(int i = 0; i < 4; ++i) {
                deg = (i*90.0 + RGSS_ANGLE) * DEG_TO_RAD;
                anti_aliased = anti_aliased + genColor(genRay(x + cos(deg)*AA_BOUND, y + sin(deg)*AA_BOUND))/4;
            }

            plot_pixel(x, y, anti_aliased.r*255.0, anti_aliased.g*255.0, anti_aliased.b*255.0);
        }
    }
    printf("Done!\n");

    fflush(stdout);
}

/*
 * Provided functions (modified for server use)
 */

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b) {
    buffer[y][x][0] = r;
    buffer[y][x][1] = g;
    buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b) {
    // Only plot jpeg for server compatibility
    plot_pixel_jpeg(x, y, r, g, b);
}

void save_jpg() {
    printf("Saving JPEG file: %s\n", filename);

    ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
    if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
        printf("Error in Saving\n");
    else
        printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found) {
    if (strcasecmp(expected, found)) {
        printf("Expected '%s ' found '%s '\n", expected, found);
        printf("Parse error, abnormal abortion\n");
        exit(0);
    }
}

void parse_doubles(FILE *file, const char *check, double p[3]) {
    char str[100];
    fscanf(file, "%s", str);
    parse_check(check, str);
    fscanf(file, "%lf %lf %lf", &p[0], &p[1], &p[2]);
    printf("%s %lf %lf %lf\n", check, p[0], p[1], p[2]);
}

void parse_rad(FILE *file, double *r) {
    char str[100];
    fscanf(file, "%s", str);
    parse_check("rad:", str);
    fscanf(file, "%lf", r);
    printf("rad: %f\n", *r);
}

void parse_shi(FILE *file, double *shi) {
    char s[100];
    fscanf(file, "%s", s);
    parse_check("shi:", s);
    fscanf(file, "%lf", shi);
    printf("shi: %f\n", *shi);
}

int loadScene(char *argv) {
    FILE *file = fopen(argv, "r");
    int number_of_objects;
    char type[50];
    Triangle t;
    Sphere s;
    Light l;
    fscanf(file, "%i", &number_of_objects);

    printf("number of objects: %i\n", number_of_objects);

    parse_doubles(file, "amb:", ambient_light);

    for (int i = 0; i < number_of_objects; i++) {
        fscanf(file, "%s\n", type);
        printf("%s\n", type);
        if (strcasecmp(type, "triangle") == 0) {
            printf("found triangle\n");
            for (int j = 0; j < 3; j++) {
                parse_doubles(file, "pos:", t.v[j].position);
                parse_doubles(file, "nor:", t.v[j].normal);
                parse_doubles(file, "dif:", t.v[j].color_diffuse);
                parse_doubles(file, "spe:", t.v[j].color_specular);
                parse_shi(file, &t.v[j].shininess);
            }

            if (num_triangles == MAX_TRIANGLES) {
                printf("too many triangles, you should increase MAX_TRIANGLES!\n");
                exit(0);
            }
            triangles[num_triangles++] = t;
        } else if (strcasecmp(type, "sphere") == 0) {
            printf("found sphere\n");

            parse_doubles(file, "pos:", s.position);
            parse_rad(file, &s.radius);
            parse_doubles(file, "dif:", s.color_diffuse);
            parse_doubles(file, "spe:", s.color_specular);
            parse_shi(file, &s.shininess);

            if (num_spheres == MAX_SPHERES) {
                printf("too many spheres, you should increase MAX_SPHERES!\n");
                exit(0);
            }
            spheres[num_spheres++] = s;
        } else if (strcasecmp(type, "light") == 0) {
            printf("found light\n");
            parse_doubles(file, "pos:", l.position);
            parse_doubles(file, "col:", l.color);

            if (num_lights == MAX_LIGHTS) {
                printf("too many lights, you should increase MAX_LIGHTS!\n");
                exit(0);
            }
            lights[num_lights++] = l;
        } else {
            printf("unknown type in scene description:\n%s\n", type);
            exit(0);
        }
    }
    return 0;
}

/*
 * Main function, modify if extra args are needed
 */

int main(int argc, char **argv) {
    if (argc != 3) {
        printf("Usage: %s <input scenefile> <output jpegname>\n", argv[0]);
        exit(0);
    }

    // Output image file
    filename = argv[2];

    // Load in abstract data to be rendered
    loadScene(argv[1]);

    // Draw the scene to global array
    draw_scene();
    // Save as a jpeg
    save_jpg();
}


#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "lib/stb/stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "lib/stb/stb_image.h"

#include "model.h"
#include "geometry.h"

int envmap_width, envmap_height;
std::vector<Vec3f> envmap;
Model duck("../duck.obj");

struct Light {
    Light(const Vec3f &p, const float i) : position(p), intensity(i) {}
    Vec3f position;
    float intensity;
};

struct Material {
    Material(const float r, const Vec4f &a, const Vec3f &color, const float spec) : refractive_index(r), albedo(a), diffuse_color(color), specular_exponent(spec) {}
    Material() : refractive_index(1), albedo(1,0,0,0), diffuse_color(), specular_exponent() {}
    float refractive_index;
    Vec4f albedo;
    Vec3f diffuse_color;
    float specular_exponent;
};

struct Sphere {
    Vec3f center;
    float radius;
    Material material;

    Sphere(const Vec3f &c, const float r, const Material &m) : center(c), radius(r), material(m) {}

    bool ray_intersect(const Vec3f &orig, const Vec3f &dir, float &t0) const {
        Vec3f L = center - orig;
        float tca = L*dir;
        float d2 = L*L - tca*tca;
        if (d2 > radius*radius) return false;
        float thc = sqrtf(radius*radius - d2);
        t0       = tca - thc;
        float t1 = tca + thc;
        if (t0 < 0) t0 = t1;
        if (t0 < 0) return false;
        return true;
    }
};

Vec3f reflect(const Vec3f &I, const Vec3f &N) {
    return I - N*2.f*(I*N);
}

Vec3f refract(const Vec3f &I, const Vec3f &N, const float eta_t, const float eta_i=1.f) { // Snell's law
    float cosi = - std::max(-1.f, std::min(1.f, I*N));
    if (cosi<0) return refract(I, -N, eta_i, eta_t); // if the ray comes from the inside the object, swap the air and the media
    float eta = eta_i / eta_t;
    float k = 1 - eta*eta*(1 - cosi*cosi);
    return k<0 ? Vec3f(1,0,0) : I*eta + N*(eta*cosi - sqrtf(k)); // k<0 = total reflection, no ray to refract. I refract it anyways, this has no physical meaning
}

bool scene_intersect(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres, Vec3f &hit, Vec3f &N, Material &material) {
    float spheres_dist = std::numeric_limits<float>::max();
    for (size_t i=0; i < spheres.size(); i++) {
        float dist_i;
        if (spheres[i].ray_intersect(orig, dir, dist_i) && dist_i < spheres_dist) {
            spheres_dist = dist_i;
            hit = orig + dir*dist_i;
            N = (hit - spheres[i].center).normalize();
            material = spheres[i].material;
        }
    }
    
    float checkerboard_dist = std::numeric_limits<float>::max();
    if (fabs(dir.y)>1e-3)  {
        float d = -(orig.y+4)/dir.y; // the checkerboard plane has equation y = -4
        Vec3f pt = orig + dir*d;
        if (d>0 && fabs(pt.x)<10 && pt.z<-10 && pt.z>-30 && d<spheres_dist) {
            checkerboard_dist = d;
            hit = pt;
            N = Vec3f(0,1,0);
            material.diffuse_color = (int(.5*hit.x+1000) + int(.5*hit.z)) & 1 ? Vec3f(.3, .3, .3) : Vec3f(.3, .2, .1);
        }
    }

    //Duck is here !
    Material      glass(1.5, Vec4f(0.0,  0.5, 0.1, 0.8), Vec3f(0.6, 0.7, 0.8),  125.);
    float duck_dist = std::numeric_limits<float>::max();
    for (size_t i=0; i < duck.nfaces(); i++) {
        float dist_i;
        //Duck before checkerboard and spheres 
        if (duck.ray_triangle_intersect(i,orig, dir, dist_i) && dist_i < duck_dist && dist_i < checkerboard_dist && dist_i < spheres_dist) {
            
            Vec3f p0 = duck.point(duck.vert(i,0));
            Vec3f p1 = duck.point(duck.vert(i,1));
            Vec3f p2 = duck.point(duck.vert(i,2));
            
            Vec3f crossed = cross((p1-p0),(p2-p0)).normalize();
            
            duck_dist = dist_i;
            hit = orig + dir*dist_i;
            N = crossed;
            material = glass;
        }
    }
    //Adding the duck_dist to the spheres and the checkboard 
    return std::min(spheres_dist, std::min(checkerboard_dist,duck_dist))<1000;
}

Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres, const std::vector<Light> &lights, size_t depth=0) {
    Vec3f point, N;
    Material material;

    if (depth>4 || !scene_intersect(orig, dir, spheres, point, N, material)) {
	double phi = atan2(dir.z, dir.x);
	double teta = acos(dir.y);

        int a = (((phi)/(2*M_PI) + .5)*envmap_width);
        int b = (teta/M_PI*envmap_height);
        return envmap[a+b*envmap_width]; // background color
    }

    Vec3f reflect_dir = reflect(dir, N).normalize();
    Vec3f refract_dir = refract(dir, N, material.refractive_index).normalize();
    Vec3f reflect_orig = reflect_dir*N < 0 ? point - N*1e-3 : point + N*1e-3; // offset the original point to avoid occlusion by the object itself
    Vec3f refract_orig = refract_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;
    Vec3f reflect_color = cast_ray(reflect_orig, reflect_dir, spheres, lights, depth + 1);
    Vec3f refract_color = cast_ray(refract_orig, refract_dir, spheres, lights, depth + 1);

    float diffuse_light_intensity = 0, specular_light_intensity = 0;
    for (size_t i=0; i<lights.size(); i++) {
        Vec3f light_dir      = (lights[i].position - point).normalize();
        float light_distance = (lights[i].position - point).norm();

        Vec3f shadow_orig = light_dir*N < 0 ? point - N*1e-3 : point + N*1e-3; // checking if the point lies in the shadow of the lights[i]
        Vec3f shadow_pt, shadow_N;
        Material tmpmaterial;
        if (scene_intersect(shadow_orig, light_dir, spheres, shadow_pt, shadow_N, tmpmaterial) && (shadow_pt-shadow_orig).norm() < light_distance)
            continue;

        diffuse_light_intensity  += lights[i].intensity * std::max(0.f, light_dir*N);
        specular_light_intensity += powf(std::max(0.f, -reflect(-light_dir, N)*dir), material.specular_exponent)*lights[i].intensity;
    }
    return material.diffuse_color * diffuse_light_intensity * material.albedo[0] + Vec3f(1., 1., 1.)*specular_light_intensity * material.albedo[1] + reflect_color*material.albedo[2] + refract_color*material.albedo[3];
}

void render(const std::vector<Sphere> &spheres, const std::vector<Light> &lights) {
	//Constantes récuperées de github
    const int eyesep = 0.2;
    const int delta = 60;
    const int   width    = 1024 + delta;
    const int   height   = 768;
    const float fov      = M_PI/3.;
    //On crée 2 buffer pour le bleu et le rouge
    std::vector<Vec3f> framebuffer1(width*height);
    std::vector<Vec3f> framebuffer2(width*height);

    #pragma omp parallel for
    for (size_t j = 0; j<height; j++) { // actual rendering loop
        for (size_t i = 0; i<width; i++) {
            float dir_x =  (i + 0.5) -  width/2.;
            float dir_y = -(j + 0.5) + height/2.;    // this flips the image at the same time
            float dir_z = -height/(2.*tan(fov/2.));
	    //We fill the buffers with an offset of eyesep 
            framebuffer1[i+j*width] = cast_ray(Vec3f(-eyesep/2,0,0), Vec3f(dir_x, dir_y, dir_z).normalize(), spheres, lights);
            framebuffer2[i+j*width] = cast_ray(Vec3f(eyesep/2,0,0), Vec3f(dir_x, dir_y, dir_z).normalize(), spheres, lights);
        }
    }

    std::vector<unsigned char> pixmap((width-delta)*height*3);
    for (size_t j = 0; j < height; ++j) {
	for (size_t i = 0; i<width-delta; i++){ 
		//On récupère les pixels qui nous intéressent
		Vec3f c1 = framebuffer1[i+delta+j*width];
		Vec3f c2 = framebuffer2[i+j*width];

		//On met l'image en niveau de gris
		float max1 = std::max(c1[0], std::max(c1[1],c1[2]));
		if (max1>1) c1 = c1*(1./max1);
            	float max2 = std::max(c2[0], std::max(c2[1], c2[2]));
            	if (max2>1) c2 = c2*(1./max2);
            	float avg1 = (c1.x+c1.y+c1.z)/3.;
            	float avg2 = (c2.x+c2.y+c2.z)/3.;

		//On remplit l'image avec les niveaux de bleu ou de rouge
            	pixmap[(j*(width-delta) + i)*3  ] = 255*avg1;
            	pixmap[(j*(width-delta) + i)*3+1] = 0;
            	pixmap[(j*(width-delta) + i)*3+2] = 255*avg2;
	} 
    }
    stbi_write_jpg("yo.jpg", width-delta, height, 3, pixmap.data(), 100);
}

int main() {
    int n = -1;
    unsigned char *pixmap = stbi_load("../envmap.jpg", &envmap_width, &envmap_height, &n, 0);
    if (!pixmap || 3!=n) {
        std::cerr << "Error: can not load the environment map" << std::endl;
        return -1;
    }
    envmap = std::vector<Vec3f>(envmap_width*envmap_height);
    for (int j = envmap_height-1; j>=0 ; j--) {
        for (int i = 0; i<envmap_width; i++) {
            envmap[i+j*envmap_width] = Vec3f(pixmap[(i+j*envmap_width)*3+0], pixmap[(i+j*envmap_width)*3+1], pixmap[(i+j*envmap_width)*3+2])*(1/255.);
        }
    }
    stbi_image_free(pixmap);

    Material      ivory(1.0, Vec4f(0.6,  0.3, 0.1, 0.0), Vec3f(0.4, 0.4, 0.3),   50.);
    Material      glass(1.5, Vec4f(0.0,  0.5, 0.1, 0.8), Vec3f(0.6, 0.7, 0.8),  125.);
    Material red_rubber(1.0, Vec4f(0.9,  0.1, 0.0, 0.0), Vec3f(0.3, 0.1, 0.1),   10.);
    Material     mirror(1.0, Vec4f(0.0, 10.0, 0.8, 0.0), Vec3f(1.0, 1.0, 1.0), 1425.);
    
    //Create a vector for the triangles (duck faces), just like for the spheres and lights.
    std::vector<Vec3i> faces;
    for(int j = 0; j < duck.nfaces(); j++){
        faces.push_back(duck.getFace(j));
    }
    
    std::vector<Sphere> spheres;
    spheres.push_back(Sphere(Vec3f(-3,    0,   -16), 2,      ivory));
    spheres.push_back(Sphere(Vec3f(-1.0, -1.5, -12), 2,      glass));
    spheres.push_back(Sphere(Vec3f( 1.5, -0.5, -18), 3, red_rubber));
    spheres.push_back(Sphere(Vec3f( 7,    5,   -18), 4,     mirror));

    std::vector<Light>  lights;
    lights.push_back(Light(Vec3f(-20, 20,  20), 1.5));
    lights.push_back(Light(Vec3f( 30, 50, -25), 1.8));
    lights.push_back(Light(Vec3f( 30, 20,  30), 1.7));

    render(spheres, lights);

    return 0;
}


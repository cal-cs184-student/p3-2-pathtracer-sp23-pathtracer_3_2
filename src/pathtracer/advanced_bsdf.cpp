#include "bsdf.h"

#include <algorithm>
#include <iostream>
#include <utility>

#include "application/visual_debugger.h"

using std::max;
using std::min;
using std::swap;

namespace CGL {

// Mirror BSDF //

    Vector3D MirrorBSDF::f(const Vector3D wo, const Vector3D wi) {
        return Vector3D();
    }

    Vector3D MirrorBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

        // TODO Project 3-2: Part 1
        // Implement MirrorBSDF
        reflect(wo, wi);
        *pdf = 1;
        return reflectance / (abs_cos_theta(*wi));
    }

    void MirrorBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Mirror BSDF"))
        {
            DragDouble3("Reflectance", &reflectance[0], 0.005);
            ImGui::TreePop();
        }
    }

// Microfacet BSDF //

    double MicrofacetBSDF::G(const Vector3D wo, const Vector3D wi) {
        return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
    }

    double MicrofacetBSDF::D(const Vector3D h) {
        // TODO Project 3-2: Part 2
        // Compute Beckmann normal distribution function (NDF) here.
        // You will need the roughness alpha.
        double e(exp(1.0));
        Vector3D n(0.0, 0.0, 1.0);
        double costheta(h.z);
        double tantheta((cross(h, n).norm()) / costheta);
        return (std::pow(e, (-(std::pow(tantheta, 2)/std::pow(alpha, 2))))/(PI*std::pow(alpha, 2)*std::pow(costheta, 4)));
    }

    Vector3D MicrofacetBSDF::F(const Vector3D wi) {
        // TODO Project 3-2: Part 2
        // Compute Fresnel term for reflection on dielectric-conductor interface.
        // You will need both eta and etaK, both of which are Vector3D.
        double costheta(wi.z);
        Vector3D etaksqr(eta * eta + k * k);
        Vector3D etatheta(2 * costheta * eta);
        Vector3D thetasqr(costheta * costheta);
        Vector3D Rs((etaksqr - etatheta + thetasqr) / (etaksqr + etatheta + thetasqr));
        Vector3D Rp((etaksqr * thetasqr - etatheta + 1) / (etaksqr * thetasqr + etatheta + 1));
        return (Rs + Rp)/2;
    }

    Vector3D MicrofacetBSDF::f(const Vector3D wo, const Vector3D wi) {
        // TODO Project 3-2: Part 2
        // Implement microfacet model here.

        Vector3D h((wi + wo) / ((wi + wo).norm()));
        Vector3D n(0.0, 0.0, 1.0);
        if ((dot(n, wi) <= 0) || (dot(n, wo) <= 0)) return Vector3D(0.0);
        return (F(wi) * G(wo, wi) * D(h))/(4.0 * dot(n, wo) * dot(n, wi));
    }

    Vector3D MicrofacetBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
        // TODO Project 3-2: Part 2
        // *Importance* sample Beckmann normal distribution function (NDF) here.
        // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
        //       and return the sampled BRDF value.
        Vector2D r(sampler.get_sample());
        double theta(atan(-std::pow(alpha, 2) * std::log(1.0 - r[0]))), phi(2 * PI * r[1]);
        Vector3D h(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
        *wi = 2.0 * (wo * h) * h - wo;
        Vector3D L_out(MicrofacetBSDF::f(wo, *wi));
        if (L_out.norm() <= EPS_D) {
            *pdf = 0;
            return L_out;
        }
        double e(exp(1.0));
        Vector3D n(0.0, 0.0, 1.0);
        double costheta(h.z);
        double sintheta(cross(h, n).norm());
        double tantheta(sintheta / costheta);
        double p_theta(((2*sintheta)/(std::pow(alpha, 2) * std::pow(costheta, 3))) * std::pow(e, (-(std::pow(tantheta, 2) / std::pow(alpha, 2)))));
        double p_phi(1 / (2 * PI));

        double p_h(p_theta * p_phi / sintheta);

        *pdf = p_h / (4 * dot(*wi, h));
        return L_out;
    }

    void MicrofacetBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Micofacet BSDF"))
        {
            DragDouble3("eta", &eta[0], 0.005);
            DragDouble3("K", &k[0], 0.005);
            DragDouble("alpha", &alpha, 0.005);
            ImGui::TreePop();
        }
    }

// Refraction BSDF //

    Vector3D RefractionBSDF::f(const Vector3D wo, const Vector3D wi) {
        return Vector3D();
    }

    Vector3D RefractionBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
        // TODO Project 3-2: Part 1
        // Implement RefractionBSDF
        *pdf = 1;
        if (!refract(wo, wi, ior)) return Vector3D();
        double eta(0.0);
        eta = ((wo.z > 0) ? (1.0 / ior): (ior));
        return transmittance / (abs_cos_theta(*wi)) / std::pow(eta, 2);
    }

    void RefractionBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Refraction BSDF"))
        {
            DragDouble3("Transmittance", &transmittance[0], 0.005);
            DragDouble("ior", &ior, 0.005);
            ImGui::TreePop();
        }
    }

// Glass BSDF //

    Vector3D GlassBSDF::f(const Vector3D wo, const Vector3D wi) {
        return Vector3D();
    }

    Vector3D GlassBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

        // TODO Project 3-2: Part 1
        // Compute Fresnel coefficient and either reflect or refract based on it.

        // compute Fresnel coefficient and use it as the probability of reflection
        // - Fundamentals of Computer Graphics page 305
        
        if (!refract(wo, wi, ior)) {
            *pdf = 1;
            reflect(wo, wi);
            return reflectance / (abs_cos_theta(*wi));
        }

        double R(0.0), n1(ior), n2(1.0);
        if (wo.z < 0) std::swap(n1, n2);
        R = std::pow(((n1 - n2) / (n1 + n2)), 2);
        if (coin_flip(R)) {
            *pdf = R;
            reflect(wo, wi);
            return R * reflectance / (abs_cos_theta(*wi));
        }
        else {
            *pdf = 1.0 - R;
            double eta(0.0);
            eta = ((wo.z > 0) ? (1.0 / ior) : (ior));
            return (1.0 - R) * transmittance / (abs_cos_theta(*wi)) / std::pow(eta, 2);
        }
    }

    void GlassBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Refraction BSDF"))
        {
            DragDouble3("Reflectance", &reflectance[0], 0.005);
            DragDouble3("Transmittance", &transmittance[0], 0.005);
            DragDouble("ior", &ior, 0.005);
            ImGui::TreePop();
        }
    }

    void BSDF::reflect(const Vector3D wo, Vector3D* wi) {

        // TODO Project 3-2: Part 1
        // Implement reflection of wo about normal (0,0,1) and store result in wi.
        *wi = 2.0 * (wo * Vector3D(0.0, 0.0, 1.0)) * Vector3D(0.0, 0.0, 1.0) - wo;

    }

    bool BSDF::refract(const Vector3D wo, Vector3D* wi, double ior) {

        // TODO Project 3-2: Part 1
        // Use Snell's Law to refract wo surface and store result ray in wi.
        // Return false if refraction does not occur due to total internal reflection
        // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
        // ray entering the surface through vacuum.
        double sign(1.0), ratio(ior);
        if (wo.z > 0) {
            sign = -sign;
            ratio = 1.0 / ratio;
        }

        double cos2_wi(1 - std::pow(ratio, 2) * (1 - std::pow(wo.z, 2)));
        if (cos2_wi < 0) {
            //*wi = Vector3D(-wo[0], -wo[1], wo[2]);
            return false;
        }

        *wi = Vector3D(-wo.x * ratio, -wo.y * ratio, sign * std::sqrt(cos2_wi)).unit();

        return true;

    }

} // namespace CGL

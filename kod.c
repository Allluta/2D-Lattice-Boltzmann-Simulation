#include <stdio.h>
#include <math.h>
#include <string.h>

#define LX (160)
#define LY (64)
#define F 9
#define cs_sq (1./3.)
#define wC (4./9.)
#define wS (1./9.)
#define wL (1./36.)

enum LatVel {_C,_N,_S,_W,_E,_NE,_NW,_SW,_SE};

double cx[F] = {0, 0, 0, -1, 1, 1, -1, -1, 1};
double cy[F] = {0, 1, -1, 0, 0, 1, 1, -1, -1};
double w[F] = {wC, wS, wS, wS, wS, wL, wL, wL, wL};

double cells[LX * LY * F];
double cells_tmp[LX * LY * F];
#define ARR(x, y, f) (f + F * (x) + F * (y) * LX)
double uLBm[LX * LY];
double vLBm[LX * LY];
double rhom[LX * LY];
#define ARRM(x, y) (x + (y) * LX)

int map[LX * LY];

double tau;
double rhoIN, uLBIN, vLBIN;
 
double fEq(int i, double rho, double u, double v) {
    double cu = (cx[i] * u + cy[i] * v) / cs_sq;
    return w[i] * rho * (1 + cu + 0.5 * cu * cu - 0.5 * (u * u + v * v) / cs_sq);
}

void collideAndStream(double *c, double *ct) {
    int x, y, i, xp, xm, yp, ym;
    double Fstar[F];
    double rho, uLB, vLB;
    double fc, fe, fw, fn, fs, fnw, fne, fse, fsw;

    for (x = 0; x < LX; x++)
        for (y = 0; y < LY; y++) {
            if (!map[ARRM(x, y)]) {
                if (x == LX - 1) {
                    for (i = 0; i < F; i++)
                        c[ARR(x, y, i)] = fEq(i, rhoIN, uLBIN, vLBIN);
                }
                if (x == 0) {
                    for (i = 0; i < F; i++)
                        c[ARR(x, y, i)] = c[ARR(x + 1, y, i)];
                }

                rho = 0;
                uLB = 0;
                vLB = 0;

                for (i = 0; i < F; i++) {
                    rho += c[ARR(x, y, i)];
                    uLB += cx[i] * c[ARR(x, y, i)];
                    vLB += cy[i] * c[ARR(x, y, i)];
                }

                uLB /= rho;
                vLB /= rho;

                rhom[ARRM(x, y)] = rho;
                uLBm[ARRM(x, y)] = uLB;
                vLBm[ARRM(x, y)] = vLB;

                fc = c[ARR(x, y, _C)];
                fs = c[ARR(x, y, _S)];
                fn = c[ARR(x, y, _N)];
                fe = c[ARR(x, y, _E)];
                fw = c[ARR(x, y, _W)];
                fne = c[ARR(x, y, _NE)];
                fnw = c[ARR(x, y, _NW)];
                fsw = c[ARR(x, y, _SW)];
                fse = c[ARR(x, y, _SE)];

                double m_rho, m_e, m_eps, m_jx, m_jy, m_qx, m_qy, m_pxx, m_pxy;
                double m_rhoe, m_ee, m_epse, m_jxe, m_jye, m_qxe, m_qye, m_pxxe, m_pxye;

                m_rho = rho;
                m_jx = rho * uLB;
                m_jy = rho * vLB;
                m_e = -(4. * fc + fe + fn + fw + fs) + 2. * (fne + fnw + fse + fsw);
                m_eps = 4. * fc - 2. * (fn + fs + fe + fw) + fnw + fne + fse + fsw;
                m_qx = 2. * (fw - fe) + fne + fse - fnw - fsw;
                m_qx = 2. * (fs - fn) + fne + fnw - fsw - fse;
                m_pxx = fe + fw - fn - fs;
                m_pxy = fne + fsw - fnw - fse;
                m_rhoe = rho;
                m_ee = -2 * rho + 3. * (m_jx * m_jx + m_jy * m_jy);
                m_epse = rho - 3. * (m_jx * m_jx + m_jy * m_jy);
                m_jxe = m_jx;
                m_qxe = -m_jx;
                m_jye = m_jy;
                m_qye = -m_jy;
                m_pxxe = (m_jx * m_jx - m_jy * m_jy);
                m_pxye = m_jx * m_jy;

                double om_e, om_eps, om_q, om_nu;
                om_e = 1.;
                om_eps = 1.;
                om_q = 1.;
                om_nu = 1. / tau;

                m_e += om_e * (m_ee - m_e);
                m_eps += om_eps * (m_epse - m_eps);
                m_qx += om_q * (m_qxe - m_qx);
                m_qy += om_q * (m_qye - m_qy);
                m_pxx += om_nu * (m_pxxe - m_pxx);
                m_pxy += om_nu * (m_pxye - m_pxy);

                Fstar[_C] = 1. / 9. * (m_rho - m_e + m_eps);
                Fstar[_E] = 1. / 36. * (4. * m_rho - m_e - 2. * m_eps) + 1. / 6. * (m_jx - m_qx) + m_pxx * .25;
                Fstar[_N] = 1. / 36. * (4. * m_rho - m_e - 2. * m_eps) + 1. / 6. * (m_jy - m_qy) - m_pxx * .25;
                Fstar[_W] = 1. / 36. * (4. * m_rho - m_e - 2. * m_eps) + 1. / 6. * (-m_jx + m_qx) + m_pxx * .25;
                Fstar[_S] = 1. / 36. * (4. * m_rho - m_e - 2. * m_eps) + 1. / 6. * (-m_jy + m_qy) - m_pxx * .25;
                Fstar[_NE] = 1. / 36. * (4. * m_rho + 2. * m_e + m_eps) + 1. / 12. * (2. * m_jx + 2. * m_jy + (m_qx + m_qy)) + m_pxy * .25;
                Fstar[_NW] = 1. / 36. * (4. * m_rho + 2. * m_e + m_eps) + 1. / 12. * (-2. * m_jx + 2. * m_jy + (-m_qx + m_qy)) - m_pxy * .25;
                Fstar[_SW] = 1. / 36. * (4. * m_rho + 2. * m_e + m_eps) - 1. / 12. * (2. * m_jx + 2. * m_jy + (m_qx + m_qy)) + m_pxy * .25;
                Fstar[_SE] = 1. / 36. * (4. * m_rho + 2. * m_e + m_eps) + 1. / 12. * (2. * m_jx - 2. * m_jy + (m_qx - m_qy)) - m_pxy * .25;

            } else {
                Fstar[_C] = c[ARR(x, y, _C)];

                Fstar[_E] = c[ARR(x, y, _W)];
                Fstar[_W] = c[ARR(x, y, _E)];
                Fstar[_S] = c[ARR(x, y, _N)];
                Fstar[_N] = c[ARR(x, y, _S)];

                Fstar[_NE] = c[ARR(x, y, _SW)];
                Fstar[_NW] = c[ARR(x, y, _SE)];
                Fstar[_SW] = c[ARR(x, y, _NE)];
                Fstar[_SE] = c[ARR(x, y, _NW)];
            }

           
            xp = (x == LX - 1 ? 0 : x + 1);
            xm = (x == 0 ? LX - 1 : x - 1);
            yp = (y == LY - 1 ? 0 : y + 1);
            ym = (y == 0 ? LY - 1 : y - 1);

            ct[ARR(x, y, _C)] = Fstar[_C];
            ct[ARR(x, yp, _N)] = Fstar[_N];
            ct[ARR(x, ym, _S)] = Fstar[_S];
            ct[ARR(xm, y, _W)] = Fstar[_W];
            ct[ARR(xp, y, _E)] = Fstar[_E];

            ct[ARR(xp, yp, _NE)] = Fstar[_NE];
            ct[ARR(xm, yp, _NW)] = Fstar[_NW];
            ct[ARR(xm, ym, _SW)] = Fstar[_SW];
            ct[ARR(xp, ym, _SE)] = Fstar[_SE];
        }
}

void setIC(double *c, double rhoi, double uLBi, double vLBi) {
    int x, y, i;
    int D = 16; 
    int radius = D / 2;
    int circle_centers[2][2] = {{LX / 3, LY / 2}, {2 * LX / 3, LY / 2}}; 
    int square_center[2] = {LX / 2, LY / 2}; 

   
    for (x = 0; x < LX; x++) {
        for (y = 0; y < LY; y++) {
            for (i = 0; i < F; i++)
                c[ARR(x, y, i)] = fEq(i, rhoi, uLBi, vLBi);
            uLBm[ARRM(x, y)] = uLBi;
            vLBm[ARRM(x, y)] = vLBi;
            rhom[ARRM(x, y)] = rhoi;
            map[ARRM(x, y)] = 0;

            
            for (int j = 0; j < 2; j++) {
                int cx = circle_centers[j][0];
                int cy = circle_centers[j][1];
                if (pow(x - cx, 2) + pow(y - cy, 2) <= pow(radius, 2)) {
                    map[ARRM(x, y)] = 1;
                }
            }

            
            int sx = square_center[0];
            int sy = square_center[1];
            if (abs(x - sx) <= radius && abs(y - sy) <= radius) {
                map[ARRM(x, y)] = 1;
            }

            if (y == 0 || y == LY - 1)
                map[ARRM(x, y)] = 1; 
        }
    }
}

void dumpStateVTK(char *fname) {
    FILE *fp;
    int x, y;

    fp = fopen(fname, "w");
    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "2D-ADE data file \n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET RECTILINEAR_GRID\n");
    fprintf(fp, "DIMENSIONS %d %d %d\n", LX, LY, 1);
    fprintf(fp, "X_COORDINATES %d double\n", LX);
    for (x = 0; x < LX; x++)
        fprintf(fp, "%lf ", (double)x + 0.5); 
    fprintf(fp, "\n");
    fprintf(fp, "Y_COORDINATES %d double\n", LY);
    for (y = 0; y < LY; y++)
        fprintf(fp, "%lf ", (double)y + 0.5); 
    fprintf(fp, "\n");
    fprintf(fp, "Z_COORDINATES 1 double\n");
    fprintf(fp, "0\n");
    fprintf(fp, "POINT_DATA %d \n", LX * LY);
    fprintf(fp, "SCALARS density double 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (y = 0; y < LY; y++)
        for (x = 0; x < LX; x++)
            fprintf(fp, "%e\n", 1 - rhom[ARRM(x, y)]);
    fprintf(fp, "SCALARS map double 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (y = 0; y < LY; y++)
        for (x = 0; x < LX; x++)
            fprintf(fp, "%d\n", map[ARRM(x, y)]);
    fprintf(fp, "\n");
    fprintf(fp, "VECTORS velocity double\n");
    for (y = 0; y < LY; y++)
        for (x = 0; x < LX; x++)
            fprintf(fp, "%e %e 0.\n", uLBm[ARRM(x, y)], vLBm[ARRM(x, y)]);
    fclose(fp);
}
int main(int argc, char **argv) {
    
    long int iter = 0;
    long int ITERMAX = 1e5;

    double Ma = 0.1;
    double Re = 1e6;
    double uLB = Ma * sqrt(cs_sq);

    double nulb = uLB * LX / Re;

    double rhoI = 1;
    double uLBI = uLB / 10;
    double vLBI = 0;

    tau = nulb / cs_sq + .5;
    rhoIN = 1.;
    uLBIN = uLB;
    vLBIN = 0;

    setIC(cells, rhoI, uLBI, vLBI);
    dumpStateVTK("state0.vtk");

    printf("%lf %lf %lf\n", uLB, nulb, tau);

    char filename[256];
    int save_interval = 10; 
    do {
        if (iter % 2)
            collideAndStream(cells_tmp, cells);
        else
            collideAndStream(cells, cells_tmp);
        if (iter % save_interval == 0) {
            sprintf(filename, "state_%04d.vtk", iter);
            dumpStateVTK(filename);
        }
        iter++;
        if (!(iter % (ITERMAX / 10))) printf("#\n");
    } while (iter < ITERMAX);

    dumpStateVTK("state.vtk");
}

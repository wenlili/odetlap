#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include "approximator.h"

#define SW 12
#define OW 6
#define RW 4
#define IT 1
#define R 0.01

// Convert coordinates to an index
inline int idx(int s, int r, int c, int nr, int nc)
{
    return s*nr*nc + r*nc + c;
}

// Process a segment
int process_segment(const float target, const int ss, const int sr, const int sc, cusp::array1d<bool,HOST> &pending, const int nslcs, const int nrows, const int ncols,
                    const cusp::array1d<float,HOST> &data, const cusp::array1d<short,HOST> &nzmask, cusp::array1d<short,HOST> &mask, cusp::array1d<float,HOST> &approximation)
{
    int count = 0;
    const int sns = nslcs/SW;
    const int snr = nrows/SW;
    const int snc = ncols/SW;
    if (!pending[idx(ss,sr,sc,snr,snc)])
        return count;
    
    const int ns0 = max(ss*SW-OW, 0);
    const int ns1 = min((ss+1)*SW+OW, nslcs);
    const int nr0 = max(sr*SW-OW, 0);
    const int nr1 = min((sr+1)*SW+OW, nrows);
    const int nc0 = max(sc*SW-OW, 0);
    const int nc1 = min((sc+1)*SW+OW, ncols);
    const int nns = ns1 - ns0;
    const int nnr = nr1 - nr0;
    const int nnc = nc1 - nc0;
    const int ds = ss*SW - ns0;
    const int dr = sr*SW - nr0;
    const int dc = sc*SW - nc0;
    
    cusp::array1d<float4,HOST> known;
    for (int i = ns0; i < ns1; i++)
        for (int j = nr0; j < nr1; j++)
            for (int k = nc0; k < nc1; k++)
                if (mask[idx(i,j,k,nrows,ncols)])
                    known.push_back(make_float4(i-ns0, j-nr0, k-nc0, data[idx(i,j,k,nrows,ncols)]));
    cusp::array1d<float,HOST> x(nns*nnr*nnc);
    
    for (int it = 0; it < IT; it++) {
        if (approximate(x, nns, nnr, nnc, known, R))
            ; //std::cout << 'c';
        else
            std::cout << '-';
        float maxe = 0;
        int sw = 0;
        int rw = 0;
        int cw = 0;
        for (int i = 0; i < SW; i++)
            for (int j = 0; j < SW; j++)
                for (int k = 0; k < SW; k++)
                    if (nzmask[idx(ss*SW+i,sr*SW+j,sc*SW+k,nrows,ncols)]) {
                        approximation[idx(ss*SW+i,sr*SW+j,sc*SW+k,nrows,ncols)] = x[idx(i+ds,j+dr,k+dc,nnr,nnc)];
                        float e = abs(x[idx(i+ds,j+dr,k+dc,nnr,nnc)] - data[idx(ss*SW+i,sr*SW+j,sc*SW+k,nrows,ncols)]);
                        if (e > maxe) {
                            maxe = e;
                            sw = i;
                            rw = j;
                            cw = k;
                        }
                    }
        if (maxe > target) {
            int iw = idx(ss*SW+sw,sr*SW+rw,sc*SW+cw,nrows,ncols);
            count++;
            mask[iw] = 1;
            known.push_back(make_float4(sw+ds, rw+dr, cw+dc, data[iw]));
            std::cout << ',' << ss << ',' << sr << ',' << sc << ',' << maxe << std::endl;
            for (int si = max(ss-1,0); si <= min(ss+1,sns-1); si++)
                for (int sj = max(sr-1,0); sj <= min(sr+1,snr-1); sj++)
                    for (int sk = max(sc-1,0); sk <= min(sc+1,snc-1); sk++)
                        if (!(si==ss && sj==sr && sk==sc)) {
                            bool c1 = ss*SW+sw >= max(si*SW-OW, 0) && ss*SW+sw < min((si+1)*SW+OW, nslcs);
                            bool c2 = sr*SW+rw >= max(sj*SW-OW, 0) && sr*SW+rw < min((sj+1)*SW+OW, nrows);
                            bool c3 = sc*SW+cw >= max(sk*SW-OW, 0) && sc*SW+cw < min((sk+1)*SW+OW, ncols);
                            if (c1 && c2 && c3)
                                pending[idx(si,sj,sk,snr,snc)] = true;
                        }
        } else {
            pending[idx(ss,sr,sc,snr,snc)] = false;
            break;
        }
    }
    
    return count;
}

// Insert a serial number into a filename that ends with ".txt"
inline std::string serial_name(char *name, float sn)
{
    std::string namestr = std::string(name);
    std::ostringstream oss;
    oss << sn;
    return namestr.substr(0, namestr.size()-4) + "-" + oss.str() + ".txt";
}

// Arguments: data nzmask nslcs nrows ncols mask values approximation
int main(int argc, char **argv)
{
    cusp::detail::timer timer = cusp::detail::timer();
    timer.start();
    
    const int nslcs = atoi(argv[3]);
    const int nrows = atoi(argv[4]);
    const int ncols = atoi(argv[5]);
    assert(!(nslcs%SW) && !(nrows%SW) && !(ncols%SW));
    const int snr   = nrows/SW;
    const int snc   = ncols/SW;
    const int size  = nslcs*nrows*ncols;
    const int ssize = size/SW/SW/SW;
    
    cusp::array1d<float,HOST> data(size);
    std::ifstream ifs(argv[1]);
    for (int i = 0; i < size; i++)
        ifs >> data[i];
    ifs.close();
    cusp::array1d<short,HOST> nzmask(size);
    ifs.open(argv[2]);
    for (int i = 0; i < size; i++)
        ifs >> nzmask[i];
    ifs.close();
    float vmin = 999999;
    float vmax = -999999;
    for (int i = 0; i < size; i++)
        if (nzmask[i]) {
            vmin = min(vmin, data[i]);
            vmax = max(vmax, data[i]);
        }
    float scale = 255/(vmax-vmin);
    for (int i = 0; i < size; i++)
        if (nzmask[i])
            data[i] = round((data[i]-vmin)*scale);
    std::cout << std::setprecision(7) << "Min: " << vmin << "\nMax: " << vmax << "\nScale: " << scale << std::endl;
    
    cusp::array1d<short,HOST> mask(size, 0);
    if (1) {
        for (int i = RW/2; i < nslcs; i += RW)
            for (int j = RW/2; j < nrows; j += RW)
                for (int k = RW/2; k < ncols; k += RW)
                    if (nzmask[idx(i,j,k,nrows,ncols)])
                        mask[idx(i,j,k,nrows,ncols)] = 1;
    } else {
        ifs.open("m.txt");
        for (int i = 0; i < size; i++)
            ifs >> mask[i];
        ifs.close();
    }
    cusp::array1d<bool,HOST> pending(ssize, true);
    cusp::array1d<float,HOST> x(size, 0);
    cusp::array1d<int,HOST> id(ssize);
    for (int i = 0; i < ssize; i++)
        id[i] = i;
    
    for (float percent = 3; percent >= 0.5; percent -= 0.5) {
        int count = 1;
        while (count) {
            count = 0;
            std::random_shuffle(id.begin(), id.end());
            for (int i = 0; i < ssize; i++)
                count += process_segment(percent/100*255, id[i]/(snr*snc), id[i]%(snr*snc)/snc, id[i]%snc, pending, nslcs, nrows, ncols, data, nzmask, mask, x);
            std::cout << "Added " << count << " known points" << std::endl;
        }
        
        std::ofstream ofs(serial_name(argv[6],percent));
        ofs << std::setprecision(7);
        for (int i = 0; i < size; i++)
            ofs << mask[i] << ((i+1)%ncols ? ' ' : '\n');
        ofs.close();
        ofs.open(serial_name(argv[7],percent));
        for (int i = 0; i < size; i++)
            if (mask[i])
                ofs << data[i] << '\n';
        ofs.close();
        ofs.open(serial_name(argv[8],percent));
        for (int i = 0; i < size; i++)
            ofs << min(max(x[i],0.0),255.0)/scale+vmin << ((i+1)%ncols ? ' ' : '\n');
        ofs.close();
        
        timer.soft_stop();
        std::cout << "Total seconds: " << timer.total_seconds() << std::endl;
        
        for (int i = 0; i < ssize; i++)
            pending[i] = true;
    }
    
    return 0;
}

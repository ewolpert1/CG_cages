import sys
import hoomd
import random
import numpy as np
from hoomd import jit
from hoomd import hpmc
from hoomd.jit import patch
from hoomd.jit import external

# Define parameter from command line arguments
ang_w=float(sys.argv[1])
ang_a=float(sys.argv[2])
jeg_w=int(sys.argv[3])
jeg_a=int(sys.argv[4])

hoomd.context.initialize('--mode=cpu')

# Define the vertices of the polyhedron
oct_vertices=[
    (-0.5, 0, 0),
    (0.5, 0, 0),
    (0, -0.5, 0),
    (0, 0.5, 0),
    (0, 0, -0.5),
    (0, 0, 0.5),
];

# Set simulation parameters 
seed = random.randint(1,1e6)
N=8
J=1
init_d = 0.2
init_a = 0.2
mr = 0.5
tuner_period = 1e3
max_part_moves = [1.0, 1.0]

system = hoomd.init.create_lattice(hoomd.lattice.sc(1.5, type_name='A'), n=N)
mc = hpmc.integrate.convex_polyhedron(seed=seed)
mc.shape_param.set('A', vertices=oct_vertices);
mc.set_params(d=init_d, a=init_a, move_ratio=mr)
tuner = hpmc.util.tune(mc, tunables=['d', 'a'], target=0.2, gamma=0.3,max_val=max_part_moves,)
r_cut =1.1; sigma_tor_w=ang_w*2; sigma_tor_a=ang_a*2; sigma_w=ang_w; sigma_a=ang_a; Jeng_a=jeg_a; Jeng_w=jeg_w

gsd_filename = filename='output/phase_%1.2f_%1.2f_%1.2f_%1.2f_%d_%d_temp.gsd'%(sigma_w*2,sigma_a*2,sigma_w,sigma_a,Jeng_w,Jeng_a)
gsd = hoomd.dump.gsd(gsd_filename, group=hoomd.group.all(), period=100000, phase=0, overwrite=True)
temp=1.1 #1.42 was used in the paper

log_filename='output/phase_%1.2f_%1.2f_%1.2f_%1.2f_%d_%d_temp.log'%(sigma_w*2,sigma_a*2,sigma_w,sigma_a,Jeng_w,Jeng_a)
logger = hoomd.analyze.log(filename=log_filename,
    quantities=['hpmc_translate_acceptance',
                'hpmc_rotate_acceptance',
                'hpmc_d',
                'hpmc_a',
                'lx',
                'hpmc_overlap_count',
                'hpmc_patch_energy',
                'time',
                'temp'],
    period=int(100),
    overwrite=False)
logger.register_callback('temp', lambda test: 1.1*0.95**(temp_+1))

for temp_ in range(20):
    temp=temp*0.95
    _temp=1/temp
    gsd.dump_state(mc)
    energy="""const float r_cut = alpha_iso[0];
const float mag = dot(r_ij, r_ij);
const float rdist = sqrt(mag);
const Scalar3 unit_vectors_w[4] = {{make_scalar3(1/sqrt(3), 1/sqrt(3), 1/sqrt(3)), 
                                 make_scalar3(1/sqrt(3), -1/sqrt(3), -1/sqrt(3)),
                                 make_scalar3(-1/sqrt(3), 1/sqrt(3), -1/sqrt(3)),
                                 make_scalar3(-1/sqrt(3), -1/sqrt(3), 1/sqrt(3))}};

const Scalar3 unit_vectors_a[4] = {{make_scalar3(-1/sqrt(3), -1/sqrt(3), -1/sqrt(3)),
                                 make_scalar3(1/sqrt(3), 1/sqrt(3), -1/sqrt(3)),
                                 make_scalar3(1/sqrt(3), -1/sqrt(3), 1/sqrt(3)),
                                 make_scalar3(-1/sqrt(3), 1/sqrt(3), 1/sqrt(3))}};

const vec3<float> norm_r = r_ij/rdist;
vec3<float> rvw_i[4], rvw_j[4], rva_i[4], rva_j[4], proj_i[3], proj_j[3];
float comp_iw[4], comp_jw[4], comp_ia[4], comp_ja[4], compa_i, compa_j, compw_i, compw_j, tor_wi[6], tor_ai[3], tor_a, tor_w;

for (int i = 0; i < 4; ++i) {{
    rvw_i[i] = rotate(q_i, vec3<float>(unit_vectors_w[i]));
    rvw_j[i] = rotate(q_j, vec3<float>(unit_vectors_w[i]));
    rva_i[i] = rotate(q_i, vec3<float>(unit_vectors_a[i]));    
    rva_j[i] = rotate(q_j, vec3<float>(unit_vectors_a[i]));
    comp_iw[i] = pow(acos(dot(norm_r, rvw_i[i])), 2);
    comp_jw[i] = pow(acos(dot(-norm_r, rvw_j[i])), 2);
    comp_ia[i] = pow(acos(dot(norm_r, rva_i[i])), 2);
    comp_ja[i] = pow(acos(dot(-norm_r, rva_j[i])), 2);    
}}

compw_i = comp_iw[0];
compw_j = comp_jw[0];
compa_i = comp_ia[0];
compa_j = comp_ja[0];
int min_index_iw=0;
int min_index_jw=0;
int min_index_ia=0;
int min_index_ja=0;
for (int i = 1; i < 4; ++i) {{
    if (comp_iw[i] < compw_i) {{
        compw_i = comp_iw[i];
        min_index_iw=i;
    }}
    if (comp_jw[i] < compw_j) {{
        compw_j = comp_jw[i];
        min_index_jw=i;
    }}
    if (comp_ia[i] < compa_i) {{
        compa_i = comp_ia[i];
        min_index_ia=i;
    }}
    if (comp_ja[i] < compa_j) {{
        compa_j = comp_ja[i];
        min_index_ja=i;
    }}    
}}


vec3<float> proj=rvw_i[0] - dot(rvw_i[0], norm_r) * norm_r;
if (min_index_iw==0)
     proj=rvw_i[2] - dot(rvw_i[2], norm_r) * norm_r;
vec3<float> proj_norm=proj/sqrt(dot(proj,proj)); 


int numj=0;
for (int i = 0; i < 4; ++i) {{
    if (i != min_index_jw) {{
        proj_j[numj] = rvw_j[i] - dot(rvw_j[i], norm_r) * norm_r;
        numj=numj+1;
    }}
}}


for (int i = 0; i < 3; ++i) {{
    proj_j[i] = proj_j[i] / sqrt(dot(proj_j[i], proj_j[i]));
    float dot_product_i_j = dot(proj_norm, proj_j[i]);
    tor_wi[i] = sqrt(pow(acos(dot_product_i_j)- 1.0472, 2));
    tor_wi[3+i] = sqrt(pow(acos(dot_product_i_j)+ 1.0472, 2));
    tor_ai[i] = sqrt(pow(acos(dot_product_i_j), 2));
}}


tor_w = tor_wi[0];
for (int i = 1; i < 6; ++i) {{
    if (tor_wi[i] < tor_w)
        tor_w = tor_wi[i];
}}


tor_a = tor_ai[0];
for (int i = 1; i < 3; ++i) {{
    if (tor_ai[i] < tor_a)
        tor_a = tor_ai[i];
}}


const float sigma_a = alpha_iso[2];
const float sigma_w = alpha_iso[4];
const float sigma_lj = 0.6;
const float Jeng_a = alpha_iso[3];
const float Jeng_w = alpha_iso[5];
const float sigma_tor_a =  sigma_a * 2;
const float sigma_tor_w =  sigma_w * 2;
const float temp = alpha_iso[8];


if (sqrt(mag) < r_cut)
  return  temp * (Jeng_a * (exp(-compa_j / (2 * pow(sigma_a, 2))) *
                          exp(-compw_i / (2 * pow(sigma_a, 2))) +
                        exp(-compw_j / (2 * pow(sigma_a, 2))) * 
                          exp(-compa_i / (2 * pow(sigma_a, 2)))) *
                (pow(sigma_lj/rdist, 12) - pow(sigma_lj/rdist, 6)) * 
                exp(-pow(tor_a, 2) / (2 * pow(sigma_tor_a, 2))) +
                Jeng_w * exp(-compw_i / (2 * pow(sigma_w, 2))) * 
                exp(-compw_j / (2 * pow(sigma_w, 2)))*
                (pow(sigma_lj/rdist, 12) - pow(sigma_lj/rdist, 6))*
                exp(-pow(tor_w, 2) / (2 * pow(sigma_tor_w, 2))));
else
    return 0.0f;
""".format(r_cut,J,sigma_a,Jeng_a,sigma_w,Jeng_w,sigma_tor_w,sigma_tor_a,temp)
    patch = hoomd.jit.patch.user(mc=mc, r_cut=r_cut, array_size=9, code=energy)
    patch.alpha_iso[:]= [r_cut,J,sigma_a,Jeng_a,sigma_w,Jeng_w,sigma_tor_w,sigma_tor_a,_temp]
    N_loops = 10001
    for i in range(N_loops):
        hoomd.run(tuner_period,quiet=True)
        tuner.update()
        print(i,'/',N_loops)


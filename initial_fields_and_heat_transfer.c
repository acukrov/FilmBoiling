#include "udf.h"
#include "sg.h"
#include "sg_mphase.h"
#include "flow.h"
#include "mem.h"
#include "metric.h"

#define SMALL 1.e-20
#define DY 6.060942e-4 
#define T_W 873.15
#define T_SAT 373.15
#define V_CELL 1.0
#define H 0.0225 

DEFINE_INIT(my_init_func, d)
{
	cell_t c;
	Thread *t;

	real xc[ND_ND];
	Domain *sub_domain;
	Domain *sd1 = DOMAIN_SUB_DOMAIN(d, 0);
	Domain *sd2 = DOMAIN_SUB_DOMAIN(d, 1);
	
	real C_storeT_radial = 0;
	real C_storeT_axial_top = 0;
	real C_storeT_axial_bottom = 0;

	int max_sub_domains;
	int phase_domain_index;
	int sub_domain_index;

	max_sub_domains = 0;
	
	sub_domain_loop(sub_domain, d, sub_domain_index)
	{
		++max_sub_domains;
	}
	
	sub_domain_loop(sub_domain, d, phase_domain_index)
	{
		thread_loop_c(t, sub_domain)
		{
			if(FLUID_THREAD_P(t))
			{
				if(DOMAIN_ID(sub_domain) == 2)
				{
					begin_c_loop_all(c, t)
					{
						C_CENTROID(xc, c, t);
						C_T(c, t) = T_SAT;
						C_U(c, t) = 0;
						C_V(c, t) = 0;						
						C_VOF(c, t) = 1;
						C_K(c, t) = 0.25;
				
						if(((xc[0] >= (-V_CELL*DY - H)) && (xc[0] <= (V_CELL*DY + H))) && ((xc[1] >= H) && (xc[1] <= (H + V_CELL*DY))))
						{
							C_VOF(c, t) = 0;
							C_K(c, t) = 0;
						}
						if((((xc[0] >= (-V_CELL*DY - H)) && (xc[0] <= -H)) || ((xc[0] >= H) && (xc[0] <= (H + V_CELL*DY)))) && ((xc[1] >= 0) && (xc[1] <= H)))
						{
							C_VOF(c, t) = 0;
							C_K(c, t) = 0; 
						}
					}
					end_c_loop_all(c, t)
				}
				else if(DOMAIN_ID(sub_domain) == 3)
				{
					int i = 0;
					int j = 0;
					int k = 0;

					begin_c_loop_all(c, t)
					{
						C_CENTROID(xc, c, t);
					
						if(((xc[0] >= -H) && (xc[0] <= H)) && ((xc[1] >= H) && (xc[1] <= (H + V_CELL*DY))) && (i==0)) 
						{
							C_storeT_radial = ((T_SAT - T_W) /log((H + V_CELL * DY)/H))*log(xc[1]/H) + T_W;
							i++;
						}
						if(((xc[0] >= (-V_CELL*DY - H)) && (xc[0] <= - H)) && ((xc[1] >= 0) && (xc[1] <= H)) && (j==0))
						{
							C_storeT_axial_bottom = ((T_SAT - T_W)/(-V_CELL*DY))*(xc[0] + H) + T_W;
							j++;	
						}
						if(((xc[0] >= H) && (xc[0] <= (H + V_CELL*DY))) && ((xc[1] >= 0) && (xc[1] <= H)) && (k==0))
						{
							C_storeT_axial_top = ((T_SAT - T_W)/(V_CELL*DY))*(xc[0] - H) + T_W;
							k++;	
						}	
		
					} 
					end_c_loop_all(c, t)
	
					begin_c_loop_all(c, t)
					{
						C_CENTROID(xc, c, t);
						C_U(c, t) = 0;
						C_V(c, t) = 0;
						C_T(c, t) = T_SAT; 	
						C_VOF(c, t) = 0;					
	
						if(((xc[0] >= -H) && (xc[0] <= H)) && ((xc[1] >= H) && (xc[1] <= (H + V_CELL*DY))))
						{
							C_T(c, t) = ((T_SAT - T_W) /log((H+ V_CELL * DY)/H))*log(xc[1]/H) + T_W;
							C_VOF(c, t) = 1;
						}
						if(((xc[0] >= (-V_CELL*DY - H)) && (xc[0] <= -H)) && ((xc[1] >= 0) && (xc[1] <= H)))
						{
							C_T(c, t) = ((T_SAT - T_W)/(-V_CELL*DY))*(xc[0] + H)  + T_W;
							C_VOF(c, t) = 1;
						}
						if(((xc[0] >= H) && (xc[0] <= (H + V_CELL*DY))) && ((xc[1] >= 0) && (xc[1] <= H)))
						{
							C_T(c, t) = ((T_SAT - T_W)/(V_CELL*DY))*(xc[0] - H) + T_W;
							C_VOF(c, t) = 1;
						}
						if(((xc[0] >= (-H - V_CELL*DY)) && (xc[0] <= -H)) && ((xc[1] >= H) && (xc[1] <= (H + V_CELL*DY))))
						{
							C_T(c, t) = (C_storeT_radial + C_storeT_axial_bottom + 2.*T_SAT)/4.;
							C_VOF(c, t) = 1;
						}
						if(((xc[0] >= H) && (xc[0] <= (H + V_CELL*DY))) && ((xc[1] >= H) && (xc[1] <= (H + V_CELL*DY)))) 
						{
							C_T(c, t) = (C_storeT_radial + C_storeT_axial_top + 2.*T_SAT)/4.;
							C_VOF(c, t) = 1;
						}
					}
					end_c_loop_all(c, t)
				}
						
			}
		}
	}
			
}

DEFINE_EXCHANGE_PROPERTY(custom_heat, cell, mix_thread, i, j)
{
	Thread *thread_l, *thread_g;
	
	/* find the threads for the liquid (primary) */
	/* and gas (secondary phases) */
	thread_l = THREAD_SUB_THREAD(mix_thread, i); /* liquid phase */
	thread_g = THREAD_SUB_THREAD(mix_thread, j); /* gas phase */
	
	real k_l;
	real d_g;
	real dx, vol, val;
	real Nu;

	k_l = C_K_L(cell, thread_l);
	d_g = C_PHASE_DIAMETER(cell, thread_g);
	vol = C_VOLUME(cell, mix_thread);
	dx = pow(vol, 1./3.);
	Nu = 2. * d_g / dx;
	val = (k_l * Nu) / d_g;
	return val;	
}

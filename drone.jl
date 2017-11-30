using JuMP
using CPLEX
using SCS
using PyPlot
include("data.jl")
function drone(instanceName::String)
 	m = Model(solver=SCSSolver())
 	K=0.2
 	v=6
 	R=20
 	w=0
 	@variable(m, b)
 	epsilon = 0.000000001  #pour avoir une inegalite stricte, on ajoute/retire � b epsilon et on fait une inegalite
 	@variable(m, u) #initialisation du terrain

	donnee=generateData(10)
	depart=donnee.start
	arrivee=donnee.destination
	gene=donnee.obstacles

	#declaration des coefficients de V(s)
		@variable(m, C[1:10])

#contrainte sur le depart
	x=depart[1]
	y=depart[2]
	t=depart[3]
 	@constraint(m, C[1]+C[2]*x+C[3]*y+C[4]*t+C[5]*x^2+C[6]*y^2+C[7]*t^2+C[8]*x*y+C[9]*x*t+C[10]*y*t <= b-epsilon)

	#contrainte sur les obstacles
	
		obs_pos=Array{Float64}(2,2)
		for i=1:length(gene)
			obs=gene[i]
			if sqrt((obs[1]-x)^2+(obs[2]-y)^2)<=R^2
		 		x_obs = obs[1]
	 			y_obs = obs[2]
				@constraints(m, begin
					obs_pos[1,1]==C[1]+C[2]*x_obs+C[3]*y_obs+C[5]*x_obs^2+C[6]*y_obs^2+C[8]*x_obs*y_obs-(b+epsilon)
					obs_pos[1,2]==(C[4]+C[9]*x_obs+C[10]*y_obs)/2
					obs_pos[2,2]==C[7]
					obs_pos[1,2]==obs_pos[2,1]
				end)
				@SDconstraint(m, obs_pos>=0)
			end
	 		i+=1
	 	end

 #contrainte sur derivee de V
#declaration des coefficients de la matrice derivee
		@variable(m, N[1:10 , 1:10] , Symmetric)
		@SDconstraint(m, N <= 0)
 			#correspondance entre coefficients par analogie


	@constraint(m, N[1,1] == C[2]*w+C[3]*v +K*C[4]*u)
 	@constraint(m, 2*N[1,2] == 2*w*C[5] + C[3]*v+K*u*C[9])
 	@constraint(m, 2*N[1,3] == 2*v*C[6]+w*C[8]+K*C[10]*u)
 	@constraint(m, 2*N[1,4] == C[9]*w-C[2]*v+C[10]*v+K*(C[7]*u-C[4]))
 	@constraint(m, 2*N[1,5] + N[4,4] == -C[9]*v-C[10]*v/2-K*C[7])
 	@constraint(m, 2*(N[1,6] + N[4,5]) == v*(C[2]/(3*2)-C[10]/2))
 	@constraint(m, 2*N[1,7] + 2*N[4,6] + N[5,5] == v*(C[3]/(4*3*2)+C[9]/(3*2)))
 	@constraint(m, 2*N[1,8]+2*N[4,7]+2*N[5,6] == v*(C[10]/(4*3*2)-C[2]/(5*4*3*2)))
 	@constraint(m, 2*N[1,9] + 2*N[4,8] +2*N[5,7]+N[6,6] == -v*(C[3]/(6*5*4*3*2)+C[9]/(5*4*3*2)))
 	@constraint(m, 2*N[1,10] + 2*N[4,9] + 2*N[5,8] + 2*N[6,7] == v*(C[2]/(7*6*5*4*3*2)-C[10]/(6*5*4*3*2)))
 	@constraint(m, 2*N[4,10]+2*N[5,9]+2*N[6,8]+N[7,7] == v*C[9]/(7*6*5*4*3*2))
 	@constraint(m, 2*N[2,4] == -2*v*C[5]-K*C[9])
 	@constraint(m, 2*2*N[2,5] == -v*C[8])
 	@constraint(m, (3*2)*N[2,6] == v*C[5])
 	@constraint(m, (4*3*2)*2*N[2,7] == v*C[8])
 	@constraint(m, (5*4*3*2)*N[2,8] == -v*C[5])
 	@constraint(m, (6*5*4*3*2)*2*N[2,9] == -v*C[8])
 	@constraint(m, (7*6*5*4*3*2)*N[2,10] == v*C[5])
 	@constraint(m, 2*N[3,4] == -v*C[8]-K*C[10])
 	@constraint(m, 2*N[3,5] == -v*C[6])
 	@constraint(m, 3*2*2*N[3,6] == v*C[8])
 	@constraint(m, 4*3*2*N[3,7] == v*C[6])
 	@constraint(m, 5*4*3*2*2*N[3,8] == -v*C[8])
 	@constraint(m, 6*5*4*3*2*N[3,9] == -v*C[6])
 	@constraint(m, 7*6*5*4*3*2*2*N[3,10] == v*C[8])

 	@constraint(m, 2*N[5,10] + 2*N[6,9] + 2*N[7,8] ==0)
 	@constraint(m, 2*N[6,10] + 2*N[7,9] + N[8,8] ==0)
 	@constraint(m, 2*N[7,10] + 2*N[8,9] ==0)
 	@constraint(m, 2*N[8,10] + N[9,9] ==0)

	@constraint(m, N[2,2]==0)
	@constraint(m, N[2,3]==0)
	@constraint(m, N[3,3]==0)
	@constraint(m, N[9,10]==0)
	@constraint(m, N[10,10]==0)

status=solve(m)
println("end")

 end

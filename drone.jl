using JuMP
using CPLEX
include("data.jl")

function drone(instanceName::String)
	m = Model(solver=CplexSolver())
	K=0.2
	v=6
	R=20
	w=0
	@variable(m, b)
	epsilon = 0.000001  #pour avoir une inegalite stricte, on ajoute/retire à b epsilon et on fait une inegalite
	@variable(m, u)


#initialisation du terrain
    file=loadDataFromFile(instanceName)
    
    x=file.start[1]
    y=file.start[2]
    t=file.start[3]
    aim=file.destination
    asteroids=file.obstacles


	for i=1:length(asteroids)
		if ((asteroids[i][1]-x)^2+(asteroids[i][2]-y)^2 <= R^2)
			push!(NotSeen, obs)		#NotSeen devrait etre une liste chainee pour etre optimal
		else
			push!(Seen, obs)		#Seen est soit un tableau, soit une liste chainee
		end
		i-=1
	end




	#declaration des coefficients de V(s)
	@variable(m, C_1)	
	@variable(m, C_2)	
	@variable(m, C_3)	
	@variable(m, C_4)	
	@variable(m, C_5)	
	@variable(m, C_6)	
	@variable(m, C_7)	
	@variable(m, C_8)	
	@variable(m, C_9)	
	@variable(m, C_10)



#contrainte sur derivee de V

	#declaration des coefficients de la matrice derivee
	@variable(m, m_11)	
	@variable(m, m_12)	
	@variable(m, m_13)	
	@variable(m, m_14)	
	@variable(m, m_15)	
	@variable(m, m_16)	
	@variable(m, m_17)	
	@variable(m, m_18)	
	@variable(m, m_19)	
	@variable(m, m_1a)	
	@variable(m, m_24)	
	@variable(m, m_25)	
	@variable(m, m_26)	
	@variable(m, m_27)	
	@variable(m, m_28)	
	@variable(m, m_29)	
	@variable(m, m_2a)	
	@variable(m, m_34)	
	@variable(m, m_35)	
	@variable(m, m_36)	
	@variable(m, m_37)	
	@variable(m, m_38)	
	@variable(m, m_39)	
	@variable(m, m_3a)	
	@variable(m, m_44)	
	@variable(m, m_45)	
	@variable(m, m_46)	
	@variable(m, m_47)	
	@variable(m, m_48)	
	@variable(m, m_49)	
	@variable(m, m_4a)	
	@variable(m, m_55)	
	@variable(m, m_56)	
	@variable(m, m_57)	
	@variable(m, m_58)	
	@variable(m, m_59)	
	@variable(m, m_66)	
	@variable(m, m_67)	
	@variable(m, m_68)	
	@variable(m, m_77)	
	@variable(m, m_5a)	
	@variable(m, m_69)	
	@variable(m, m_6a)	
	@variable(m, m_78)	
	@variable(m, m_79)	
	@variable(m, m_7a)	
	@variable(m, m_88)	
	@variable(m, m_89)	
	@variable(m, m_8a)	
	@variable(m, m_99)
		
	#correspondance entre coefficients par analogie
	@constraint(m, m_11 == C_2*w+C_3*v +K*C_4*u)
	@constraint(m, 2*m_12 == 2*w*C_5 + C_3*v+K*u*C_9)
	@constraint(m, 2*m_13 == 2*v*C_6+w*C_8+K*C_10*u)
	@constraint(m, 2*m_14 == C_9*w-C_2*v+C_10*v+K*(C_7*u-C_4))
	@constraint(m, 2*m_15 + m_44 == -C_9*v-C_10*v/2-K*C_7)
	@constraint(m, 2*(m_16 + m_45) == v*(C_2/(3*2)-C_10/2))
	@constraint(m, 2*m_17 + 2*m_46 + m_55 == v*(C_3/(4*3*2)+C_9/(3*2)))
	@constraint(m, 2*m_18+2*m_47+2*m_56 == v*(C_10/(4*3*2)-C_2/(5*4*3*2)))
	@constraint(m, 2*m_19 + 2*m_48 +2*m_57+m_66 == -v*(C_3/(6*5*4*3*2)+C_9/(5*4*3*2)))
	@constraint(m, 2*m_1a + 2*m_49 + 2*m_58 + 2*m_67 == v*(C_2/(7*6*5*4*3*2)-C_10/(6*5*4*3*2)))
	@constraint(m, 2*m_4a+2*m_59+2*m_68+m_77 == v*C_9/(7*6*5*4*3*2))
	@constraint(m, 2*m_24 == -2*v*C_5-K*C_9)
	@constraint(m, 2*2*m_25 == -v*C_8)
	@constraint(m, (3*2)*m_26 == v*C_5)
	@constraint(m, (4*3*2)*2*m_27 == v*C_8)
	@constraint(m, (5*4*3*2)*m_28 == -v*C_5)
	@constraint(m, (6*5*4*3*2)*2*m_29 == -v*C_8)
	@constraint(m, (7*6*5*4*3*2)*m_2a == v*C_5)
	@constraint(m, 2*m_34 == -v*C_8-K*C_10)
	@constraint(m, 2*m_35 == -v*C_6)
	@constraint(m, 3*2*2*m_36 == v*C_8)
	@constraint(m, 4*3*2*m_37 == v*C_6)
	@constraint(m, 5*4*3*2*2*m_38 == -v*C_8)
	@constraint(m, 6*5*4*3*2*m_39 == -v*C_6)
	@constraint(m, 7*6*5*4*3*2*2*m_3a == v*C_8)
	
	@constraint(m, 2*m_5a + 2*m_69 + 2*m_78 ==0)
	@constraint(m, 2*m_6a + 2*m_79 + m_88 ==0)
	@constraint(m, 2*m_7a + 2*m_89 ==0)
	@constraint(m, 2*m_8a + m_99 ==0) 

	A = [ m_11 m_12 m_13 m_14 m_15 m_16 m_17 m_18 m_19 m_1a; 
	      m_12 0    0    m_24 m_25 m_26 m_27 m_28 m_29 m_2a;
	      m_13 0    0    m_34 m_35 m_36 m_37 m_38 m_39 m_3a;
	      m_14 m_24 m_34 m_44 m_45 m_46 m_47 m_48 m_49 m_4a;
	      m_15 m_25 m_35 m_45 m_55 m_56 m_57 m_58 m_59 m_5a;
	      m_16 m_26 m_36 m_46 m_56 m_66 m_67 m_68 m_69 m_6a;
	      m_17 m_27 m_37 m_47 m_57 m_67 m_77 m_78 m_79 m_7a;
	      m_18 m_28 m_38 m_48 m_58 m_68 m_78 m_88 m_89 m_8a;
	      m_19 m_29 m_39 m_49 m_59 m_69 m_79 m_89 m_99 0   ;
	      m_1a m_2a m_3a m_4a m_5a m_6a m_7a m_8a 0    0   ]
	@SDconstraint(m, A<=0)

#contrainte sur le depart
	@constraint(m, C_1+C_2*x+C_3*y+C_4*t+C_5*x^2+C_6*y^2+C_7*t^2+C_8*x*y+C_9*x*t+C_10*y*t <= b-epsilon)

#contrainte sur les obstacles
	
	for i=1:length(Seen)
		x_obs = Seen[i][1]
		y_obs = Seen[i][2]
		@constraint(m,C_1-b+C_2*x_obs+C_3*y_obs+C_5*x_obs^2+C_6*y_obs^2+C_8*x_obs*y_obs-(b+epsilon) >= 0)				# V(0) >= b+epsilon
		delta = (C_9*x_obs+C_10*y_obs)^2-4*C_7*(C_1-b+C_2*x_obs+C_3*y_obs+C_5*x_obs^2+C_6*y_obs^2+C_8*x_obs*y_obs-(b+epsilon))
		if delta >=0
			s1=(-(C_9*x_obs+C_10*y_obs)+sqrt(delta))/(2*C_7)							#si il y a un zero sur la fonction
			s2=(-(C_9*x_obs+C_10*y_obs)-sqrt(delta))/(2*C_7)							#on ajoute la contrainte qu'il soit en dehors
			s1=min(abs(s1), abs(s2))										#de [-pi pi]
			@constraint(m, s1 >= 3.141592)
		end
		i+=1
	end


#on a tenu compte de toutes les contraintes, on doit resoudre
	status=solve(m)
end

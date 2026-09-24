-- Robert Schaefer 23-09-2026
-- Simulation of flow over cone at Mach=1.5
--

-- global settings
config.solver_mode = 'transient'
config.dimensions = 2
config.axisymmetric = true

-- gas model and flow conditions
setGasModel('ideal-air.gas')

p_inf = 95.84e3		-- Pa
T_inf = 1103 		-- K
u_inf = 1000.0		-- m/s

inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf}
initial = FlowState:new{p=0.2*p_inf, T=T_inf, velx=0.0}

-- geometry and grid
a0 = Vector3:new{x=0.0, y=0.0}
b0 = Vector3:new{x=0.2, y=0.0}
c0 = Vector3:new{x=1.0, y=0.29118}
a1 = Vector3:new{x=0.0, y=1.0}
b1 = Vector3:new{x=0.2, y=1.0}
c1 = Vector3:new{x=1.0, y=1.0}

patch0 = CoonsPatch:new{p00=a0, p10=b0, p11=b1, p01=a1}
patch1 = AOPatch:new{p00=b0, p10=c0, p11=c1,  p01=b1}

grid0 = registerFluidGrid{
	grid=StructuredGrid:new{psurface=patch0, niv=11, njv=41},
	fsTag='inflow',
	bcTags={west='inflow'}
}
		
grid1 = registerFluidGrid{
	grid=StructuredGrid:new{psurface=patch1, niv=41, njv=41},
	fsTag='initial',
	bcTags={east='outflow'}
}

identifyGridConnections()

-- setup fluid blocks
flowDict = {
	initial=initial,
	inflow=inflow,
}
bcDict = {
	inflow=InFlowBC_Supersonic:new{flowState=inflow},
	outflow=OutFlowBC_Simple:new{}
}	

makeFluidBlocks(bcDict, flowDict)

-- simulations settings
config.max_time = 5.0e-3 --s
config.max_step = 3000
config.cfl_value = 0.5
config.dt_plot = 1.5e-3
config.dt_init = 1.0e-6






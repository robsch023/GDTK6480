-- Robert Schaefer 23-09-2026
-- Simulation of flow over sphere at Mach 6.0
--

-- global settings
config.solver_mode = 'transient'
config.dimensions = 2
config.axisymmetric = true

-- gas model and flow conditions
nsp, nmodes, gm = setGasModel('ideal-Argon.gas')

p_inf = 10e3		-- Pa
T_inf = 300			-- K
M_inf = 6.0			

-- find velocity from gas conditions and sound speed
gs_inf = GasState:new{gm}
gs_inf.p = p_inf
gs_inf.T = T_inf
gm:updateThermoFromPT(gs_inf)
gm:updateSoundSpeed(gs_inf)

u_inf = M_inf * gs_inf.a -- freestream velocity for Mach 6

inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf}
initial = FlowState:new{p=0.2*p_inf, T=T_inf, velx=0.0}

r_in  = 0.01
r_out = 0.03
centre0 = Vector3:new{x=0.0, y=0.0}
centre1 = Vector3:new{x=0.01, y=0.0}

-- inner arc points
p_in_180  = Vector3:new{x=-r_in, y=0.0}
p_in_90   = Vector3:new{x=0.0,   y=r_in}

-- outer arc points
p_out_180 = Vector3:new{x=centre1.x - r_out, y=centre1.y}
p_out_90 = Vector3:new{x=0.0, y=math.sqrt(r_out^2 - centre1.x^2)} 

south = Arc:new{p0=p_in_180,  p1=p_in_90,  centre=centre0}
north = Arc:new{p0=p_out_180, p1=p_out_90, centre=centre1}
west  = Line:new{p0=p_in_180, p1=p_out_180}
east  = Line:new{p0=p_in_90,  p1=p_out_90}

patch0 = CoonsPatch:new{north=north, east=east, south=south, west=west}

grid0 = StructuredGrid:new{psurface=patch0, niv=61, njv=41}

registerFluidGrid{
   grid = grid0,
   fsTag = "initial",
   bcTags = {north="inflow", east="outflow"}
}

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
config.max_time = 1.0e-3 --s
config.max_step = 100000
config.cfl_value = 0.5
config.dt_plot = 5.0e-5
config.dt_init = 1.0e-7

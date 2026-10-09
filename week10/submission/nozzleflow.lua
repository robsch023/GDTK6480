-- Robert Schaefer 09-10-2026
-- Simulation of flow through a nozzle
--

-- global settings
config.solver_mode = 'transient'
config.dimensions = 2
config.axisymmetric = true

setGasModel('ideal-air.gas')

-- res. and initial low-p states
p_r = 500.0		-- Pa
T_r = 300.0		-- K
reservoir = FlowState:new{p=p_r, T=T_r}
lowP      = FlowState:new{p=5.0, T=T_r}

-- geometry
wall = Spline2:new{filename="nozzle-contour.dat"}

-- corners
b = wall(0.0)   -- inlet, on the wall
c = wall(1.0)   -- exit, on the wall
a = Vector3:new{x=b.x, y=0.0}   -- inlet, on the axis
d = Vector3:new{x=c.x, y=0.0}   -- exit, on the axis

west  = Line:new{p0=a, p1=b}
east  = Line:new{p0=d, p1=c}
south = Line:new{p0=a, p1=d}
north = wall

x_throat, y_min = b.x, b.y
for i = 0, 1000 do
   local p = wall(i/1000)
   if p.y < y_min then y_min = p.y; x_throat = p.x end
end

patch = CoonsPatch:new{north=north, east=east, south=south, west=west}
nx, ny = 200, 40
grid0  = StructuredGrid:new{psurface=patch, niv=nx+1, njv=ny+1}

-- initial condition as a function of position
function initialState(x, y, z)
   if x <= x_throat then return reservoir else return lowP end
end

-- grid
registerFluidGrid{
   grid = grid0,
   fsTag = "initial",
   bcTags = {west="inflow", east="outflow", north="wall", south="wall"}
}

flowDict = { initial=initialState }
bcDict = {
   inflow  = InFlowBC_FromStagnation:new{stagnationState=reservoir},
   outflow = OutFlowBC_Simple:new{},
   wall    = WallBC_WithSlip:new{}
}
makeFluidBlocks(bcDict, flowDict)

-- simulation settings
config.max_time = 5.0e-3
config.max_step = 1000000
config.cfl_value = 0.5
config.dt_init = 1.0e-7
config.dt_plot = 4.0e-5
config.dt_history = 1.0e-5

-- pressure transducers
setHistoryPoint{x=b.x, y=b.y}
setHistoryPoint{x=c.x, y=c.y}

-- Robert Schaefer 12-10-26
-- Preliminary simulation of 2 ramp inlet
-- Testing broader setup parameters on course grid
--

-- global settings
config.solver_mode = 'transient'
config.dimensions = 2
config.axisymmetric = false

-- gas model and flow conditions
nsp, nmodes, gm = setGasModel('ideal-air.gas')

p_inf = 1249		-- Pa
T_inf = 227.0		-- K - estimate from altitude assoc. with given pressure
M_inf = 7.44			

-- find velocity from gas conditions
local gamma, R = 1.4, 287.0 -- gamma given, R=287 J/kgK for dry air
u_inf = M_inf*math.sqrt(gamma*R*T_inf) -- free stream velocity

inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf}
initial = inflow

-- GEOMETRY
theta1 = 6.16        -- deg, first ramp angle from freestream direction
theta2 = 7.61        -- deg, additional turning at second ramp

A = Vector3:new{x=0.000, y=0.361}
B = Vector3:new{x=0.994, y=0.253}
C = Vector3:new{x=1.879, y=0.037}
L_const = 0.30       -- m, length of constant-area section

-- RESOLUTIONS
rf  = 1
local function nc(n) return math.max(2, math.floor(n*rf + 0.5)) end
nx0 = nc(50)         -- block 0: A -> B
nx1 = nc(50)         -- block 1: B -> C
nx2 = nc(30)         -- block 2
ny  = nc(24)

-- GEOMETRY
O0 = Vector3:new{x=A.x,           y=0.0}      -- symmetry line below A
O1 = Vector3:new{x=B.x,           y=0.0}      -- below B
O2 = Vector3:new{x=C.x,           y=0.0}      -- below C
D  = Vector3:new{x=C.x+L_const,   y=0.0}      -- end of duct, symmetry line
E  = Vector3:new{x=C.x+L_const,   y=C.y}      -- end of duct, top wall

patch0 = CoonsPatch:new{p00=O0, p10=O1, p11=B, p01=A}
patch1 = CoonsPatch:new{p00=O1, p10=O2, p11=C, p01=B}
patch2 = CoonsPatch:new{p00=O2, p10=D,  p11=E, p01=C}

grid0 = StructuredGrid:new{psurface=patch0, niv=nx0+1, njv=ny+1}
grid1 = StructuredGrid:new{psurface=patch1, niv=nx1+1, njv=ny+1}
grid2 = StructuredGrid:new{psurface=patch2, niv=nx2+1, njv=ny+1}


-- boundary conditions and block registration (may need to add interfaces)
flowDict = {initial=initial, inflow=inflow}
bcDict = {
   inflow   = InFlowBC_Supersonic:new{flowState=inflow},
   outflow  = OutFlowBC_Simple:new{},
   wall     = WallBC_WithSlip:new{},
   symmetry = WallBC_WithSlip:new{},      
}

registerFluidGrid{grid=grid0, fsTag="initial",
                  bcTags={west="inflow", south="symmetry", north="wall"}}
registerFluidGrid{grid=grid1, fsTag="initial",
                  bcTags={south="symmetry", north="wall"}}
registerFluidGrid{grid=grid2, fsTag="initial",
                  bcTags={south="symmetry", north="wall", east="outflow"}}
identifyGridConnections()
makeFluidBlocks(bcDict, flowDict)

-- simulations settings
config.max_time = 6.0e-3 --s
config.max_step = 100000
config.cfl_value = 0.5
config.dt_plot = 5.0e-4
config.dt_init = 1.0e-7

--HISTORY POINTS for temporal convergence checks
-- wall height at position x
local function y_wall(x)
   if x <= B.x then
      return A.y + (B.y - A.y)*(x - A.x)/(B.x - A.x)
   elseif x <= C.x then
      return B.y + (C.y - B.y)*(x - B.x)/(C.x - B.x)
   else
      return C.y
   end
end
 
-- probe at fraction f of the local wall height
local function probe(x, f)
   setHistoryPoint{x=x, y=f*y_wall(x)}
end
 
-- probes behind ramp shocks
probe(0.5*(A.x + B.x), 0.8)     -- ramp 1
probe(0.5*(B.x + C.x), 0.8)     -- ramp 2
 
-- downstream from c
x_thr = C.x + 0.02
probe(x_thr, 0.5)
 
-- constant-area section
probe(C.x + 0.5*L_const, 0.5)
probe(C.x + 0.9*L_const, 0.5)

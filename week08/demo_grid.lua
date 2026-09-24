-- Author: RJG
-- Date: 2024-09-11
--
-- Week 08 example grid building
--
-- To run this file:
--
-- > e4shared --custom-script --script-file=demo_grid.lua


-- points for construction
a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=1.0, y=0.0}
c = Vector3:new{x=2.0, y=1.0}
d = Vector3:new{x=0.0, y=2.0}

ctr = Vector3:new{x=2.0, y=0.0}
ad = Vector3:new{x=0.7, y=1.0}
dc = Vector3:new{x=0.7, y=2.5}

-- paths for edges of domain
south = Line:new{p0=a, p1=b}
east = Arc:new{p0=b, p1=c, centre=ctr}
north = Bezier:new{points={d, dc, c}}
west = Bezier:new{points={a, ad, d}}

-- mathematical definition as parametric surface
patch = CoonsPatch:new{north=north, east=east, south=south, west=west}

-- grid definition (includes discretisation)
grid = StructuredGrid:new{psurface=patch, niv=21, njv=21}

grid:write_to_vtk_file('my-grid.vtk')


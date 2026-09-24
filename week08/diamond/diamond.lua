-- WEEK08 SUBMISSION
-- Grid generation - Diamond

a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=-1.0, y=0.2}
c = Vector3:new{x=1.0, y=0.2}
d = Vector3:new{x=0.0, y=1.0}

ab = Line:new{p0=a, p1=b}
ac = Line:new{p0=a, p1=c}
bd = Line:new{p0=b, p1=d}
cd = Line:new{p0=c, p1=d}

patch0 = CoonsPatch:new{north=cd, east=bd, south=ab, west=ac}

grid0 = StructuredGrid:new{psurface=patch0, niv=21, njv=11}
grid0:write_to_vtk_file('diamond.vtk')

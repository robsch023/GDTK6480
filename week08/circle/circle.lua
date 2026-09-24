-- WEEK08 SUBMISSION
-- Grid generation - circle
R = 1.0

a = Vector3:new{x=R*math.cos(math.pi/4), y=R*math.sin(math.pi/4)}
b = Vector3:new{x=R*math.cos(3*math.pi/4), y=R*math.sin(3*math.pi/4)}
c = Vector3:new{x=R*math.cos(5*math.pi/4), y=R*math.sin(5*math.pi/4)}
d = Vector3:new{x=R*math.cos(7*math.pi/4), y=R*math.sin(7*math.pi/4)}
centre = Vector3:new{x=0.0, y=0.0}

north = Arc:new{p0=b, p1=a, centre=centre}
south = Arc:new{p0=c, p1=d, centre=centre}
east = Arc:new{p0=d, p1=a, centre=centre}
west = Arc:new{p0=c, p1=b, centre=centre}

patch0 = CoonsPatch:new{north=north, east=east, south=south, west=west}
grid0 = StructuredGrid:new{psurface=patch0, niv=21, njv=21}
grid0:write_to_vtk_file('circle.vtk')

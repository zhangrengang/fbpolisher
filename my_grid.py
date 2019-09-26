import drmaa
with drmaa.Session() as s:
	name = s.drmaaImplementation
print [name]

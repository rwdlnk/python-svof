import vtk

fname = "/media/rod/ResearchIII/ResearchIII/githubRepos/svof/src/master/f77/slimMaster/python_refactor/slimMaster_test/RT160x200-100.vtk"

r = vtk.vtkDataSetReader()
r.SetFileName(fname)
r.ReadAllScalarsOn()
r.Update()

data = r.GetOutput()
pd = data.GetPointData()
print("Arrays:", pd.GetNumberOfArrays())
for i in range(pd.GetNumberOfArrays()):
    print(i, pd.GetArrayName(i))


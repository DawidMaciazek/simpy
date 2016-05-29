from simpy import analyze
from simpy import write


def xyz_test(ifile_name):
    traj = analyze.Traj(ifile_name, "xyz")
    frame = traj.read()
    print frame['coord']
    print frame['comment']

def lammpstrj_test(ifile_name):
    traj = analyze.Traj(ifile_name, "lammpstrj")
    frame = traj.read()
    print frame['velocity']

def test_read_write(ifile_name, ofile_name):
    traj = analyze.Traj(ifile_name, "lammpstrj")
    frame = traj.read()
    wxyz = write.xyz(ofile_name, 0)
    wxyz.write(frame)






#xyz_test("../resources/surf_test.xyz")
test_read_write('../resources/cluster.lammpstrj', 'out.xyz')

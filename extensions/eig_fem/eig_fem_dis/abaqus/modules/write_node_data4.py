import numpy as np

def write_node_data4(fname2, rn, links, xyzlimit):
    fnamedata = fname2 + '.data'

    version_number = 5
    filesegments_number = 1
    domainID = 0

    xBoundMin = -xyzlimit
    yBoundMin = -xyzlimit
    zBoundMin = -xyzlimit
    xBoundMax = xyzlimit
    yBoundMax = xyzlimit
    zBoundMax = xyzlimit

    Nnodes = np.sum(rn[:, 3] != -1)

    with open(fnamedata, 'w') as fid:
        fid.write(f"dataFileVersion =   {version_number}\n")
        fid.write(f"numFileSegments =   {filesegments_number}\n")
        fid.write("minCoordinates = [\n  {:.6e}\n  {:.6e}\n  {:.6e}\n  ]\n".format(xBoundMin, yBoundMin, zBoundMin))
        fid.write("maxCoordinates = [\n  {:.6e}\n  {:.6e}\n  {:.6e}\n  ]\n".format(xBoundMax, yBoundMax, zBoundMax))
        fid.write(f"nodeCount =   {Nnodes}\n")
        fid.write("dataDecompType =   2\n")
        fid.write("dataDecompGeometry = [\n  1\n  1\n  1\n  ]\n")
        fid.write("\n#\n#  END OF DATA FILE PARAMETERS\n#\n\n")
        fid.write("domainDecomposition = \n")
        fid.write("# Dom_ID  Minimum XYZ bounds   Maximum XYZ bounds\n")
        fid.write("  {}  {:12.6e} {:12.6e} {:12.6e}    {:12.6e}   {:12.6e}  {:12.6e}\n".format(
            domainID, xBoundMin, yBoundMin, zBoundMin, zBoundMax, yBoundMax, xBoundMax))

        fid.write("nodalData = \n")
        fid.write("#  Primary lines: node_tag, x, y, z, num_arms, constraint\n")
        fid.write("#  Secondary lines: arm_tag, burgx, burgy, burgz, nx, ny, nz\n")

        neighbors = np.zeros(rn.shape[0], dtype=int)
        for i in range(rn.shape[0]):
            neighbors[i] = np.sum(links[:, 0] == i + 1) + np.sum(links[:, 1] == i + 1)
            if neighbors[i] != int(rn[i, 5]):
                print(f"    + Arm number {int(rn[i,5])} in node {int(rn[i,4])} is incorrect. ({neighbors[i]} is correct!)")
                rn[i, 5] = neighbors[i]

        for i in range(rn.shape[0]):
            if rn[i, 5] == 1 and rn[i, 3] == 0:
                print(f"    + Node flag {int(rn[i,3])} in node {int(rn[i,4])} is incorrect.")
                rn[i, 3] = 6

        for i in range(rn.shape[0]):
            fid.write(" {},{} {:12.10g} {:12.10g} {:12.10g} {:d} {:d}\n".format(
                domainID, i, rn[i, 0], rn[i, 1], rn[i, 2], int(rn[i, 5]), int(rn[i, 3])))

            for j in range(links.shape[0]):
                if int(links[j, 0]) == i + 1:
                    fid.write("   {},{} {:12.10e} {:12.10e} {:12.10e}\n".format(
                        domainID, int(links[j, 1]) - 1, links[j, 2], links[j, 3], links[j, 4]))
                    fid.write("       {:12.10e} {:12.10e} {:12.10e}\n".format(
                        links[j, 5], links[j, 6], links[j, 7]))
            for j in range(links.shape[0]):
                if int(links[j, 1]) == i + 1:
                    fid.write("   {},{} {:12.10e} {:12.10e} {:12.10e}\n".format(
                        domainID, int(links[j, 0]) - 1, -links[j, 2], -links[j, 3], -links[j, 4]))
                    fid.write("       {:12.10e} {:12.10e} {:12.10e}\n".format(
                        links[j, 5], links[j, 6], links[j, 7]))

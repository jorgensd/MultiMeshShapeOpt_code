import sys
if __name__=="__main__":
    try:
        res = float(sys.argv[1])
    except:
        print("Invalid resolution")
        sys.exit(1)
    with open('meshes/cable.geo', 'r') as file:
        data = file.readlines()
        data[1] = "res = %s;\n" % res
        
    with open('meshes/cable.geo', 'w') as file:
        file.writelines( data )

    with open('meshes/inner_cable_halo.geo', 'r') as file:
        data = file.readlines()
        data[1] = "res = %s;\n" % res
        
    with open('meshes/inner_cable_halo.geo', 'w') as file:
        file.writelines( data )
 
    import os
    # os.system('trash meshes/cable.msh')
    os.system('make -C meshes')

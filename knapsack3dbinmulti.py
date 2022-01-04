# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 10:39:47 2021

@author: harres.tariq
"""
from ortools.linear_solver import pywraplp
import numpy as np
import matplotlib.pyplot as plt
# %%
def create_data_model():
    """Create the data for the example."""
    data = {}
    widths = [3,1,2,1,2,2]
    lengths = [3,1,1,3,2,2]
    heights=[1,1,2,2,2,2]
    values = [widths[i]*lengths[i]*heights[i] for i in range(len(widths))]#[1 for i in range(len(widths))]
    data['widths'] = widths
    data['lengths']=lengths
    data['heights']=heights
    data['values'] = values
    data['items'] = list(range(len(widths)))
    data['num_items'] = len(widths)
    num_bins = 2
    data['bins'] = list(range(num_bins))
    data['bin_max_width'] = [3,3]
    data['bin_max_length'] = [3,3]
    data['bin_max_height']=[3,3]
    return data
# %%
def make_boxes(data,cube_origins,binj):
    x,y,z=np.indices((data['bin_max_width'][binj],data['bin_max_length'][binj],data['bin_max_height'][binj]))
    voxels=(np.zeros((data['bin_max_width'][binj],data['bin_max_length'][binj],data['bin_max_height'][binj]))==True)
    colors = np.empty(voxels.shape, dtype=object)
    #n=len(cube_origins)
    #clr = cm.rainbow(np.linspace(0, 1, n))
    clr=['firebrick','yellow','green', 'red','blue','magenta','pink','cyan','khaki','gold','rosybrown','lightcoral','tan','olive']
       
    for i in data['items']:
        try:
            if(cube_origins[binj]['x'][i]>=0):
                cube=(x >= cube_origins[binj]['x'][i]) & (x <= cube_origins[binj]['x'][i]+data['widths'][i]-1) & (y >= cube_origins[binj]['y'][i]) & (y <= cube_origins[binj]['y'][i]+data['lengths'][i]-1) & (z >= cube_origins[binj]['z'][i]) & (z <= cube_origins[binj]['z'][i]+data['heights'][i]-1)
                colors[cube]=clr[i]
                voxels |= cube
        except:
            continue

    ax = plt.figure().add_subplot(projection='3d')
    ax.voxels(voxels, facecolors=colors, edgecolor='k')

    plt.show()
# %%
def main():
    data = create_data_model()

    # Create the mip solver with the SCIP backend.
    solver = pywraplp.Solver.CreateSolver('SCIP')

    x={}
    y={}    
    z={}
    ind = {}
    for i in data['items']:
        for j in data['bins']:
            x[(i,j)]=solver.IntVar(0, data['bin_max_width'][j],'x'+str(i)+str(j))
            y[(i,j)]=solver.IntVar(0, data['bin_max_length'][j],'y'+str(i)+str(j))
            z[(i,j)]=solver.IntVar(0, data['bin_max_height'][j],'z'+str(i)+str(j))
            ind[(i, j)] = solver.IntVar(0, 1, 'ind_%i_%i' % (i, j))
    
    x_bool={}
    y_bool={}
    z_bool={}
    M=int(1e6)
    for i in data['items']:
        for j in data['bins']:
            for k in data['items']:
                if i>=k:
                    continue
                else:
                    x_bool[(i,j,k)] = solver.IntVar(0,1,"x_bool_"+str(i)+str(j)+str(k))
                    y_bool[(i,j,k)] = solver.IntVar(0,1,"y_bool_"+str(i)+str(j)+str(k))
                    z_bool[(i,j,k)] = solver.IntVar(0,1,"z_bool_"+str(i)+str(j)+str(k))
                    x_bool[(k,j,i)] = solver.IntVar(0,1,"x_bool_"+str(k)+str(j)+str(i))
                    y_bool[(k,j,i)] = solver.IntVar(0,1,"y_bool_"+str(k)+str(j)+str(i))
                    z_bool[(k,j,i)] = solver.IntVar(0,1,"z_bool_"+str(k)+str(j)+str(i))
                    solver.Add(x[(i,j)]+data['widths'][i]<=x[(k,j)]+M*x_bool[(i,j,k)])
                    solver.Add(y[(i,j)]+data['lengths'][i]<=y[(k,j)]+M*y_bool[(i,j,k)])
                    solver.Add(z[(i,j)]+data['heights'][i]<=z[(k,j)]+M*z_bool[(i,j,k)])
                    solver.Add(x[(k,j)]+data['widths'][k]<=x[(i,j)]+M*x_bool[(k,j,i)])
                    solver.Add(y[(k,j)]+data['lengths'][k]<=y[(i,j)]+M*y_bool[(k,j,i)])
                    solver.Add(z[(k,j)]+data['heights'][k]<=z[(i,j)]+M*z_bool[(k,j,i)])
                    tmp=[]
                    tmp.append(x_bool[(i,j,k)])
                    tmp.append(y_bool[(i,j,k)])
                    tmp.append(z_bool[(i,j,k)])
                    tmp.append(x_bool[(k,j,i)])
                    tmp.append(y_bool[(k,j,i)])
                    tmp.append(z_bool[(k,j,i)])
                    solver.Add(sum(tmp)<=5)

        for j in data['bins']:
            solver.Add(x[(i,j)]+data['widths'][i]<=data['bin_max_width'][j]+M*(1-ind[(i,j)]))
            solver.Add(y[(i,j)]+data['lengths'][i]<=data['bin_max_length'][j]+M*(1-ind[(i,j)]))
            solver.Add(z[(i,j)]+data['heights'][i]<=data['bin_max_height'][j]+M*(1-ind[(i,j)]))
    
    for i in data['items']:
        solver.Add(sum(ind[(i, j)] for j in data['bins']) <= 1)
        
    for j in data['bins']:    
        solver.Add(sum(ind[(i,j)]*data['values'][i] for i in data['items'])<=
                   data['bin_max_width'][j]*data['bin_max_length'][j]*data['bin_max_height'][j])
    

    # Objective
    objective = solver.Objective()

#    for i in data['items']:
#        for j in data['bins']:
#            objective.SetCoefficient(ind[(i,j)], data['values'][i])
#    objective.SetMaximization()
    C=0.01
#    D=0.001
#    E=0.0001
    objective_terms = []
    for i in data['items']:
        for j in data['bins']:
            objective_terms.append((len(data['bins'])-j)*data['values'][i]*ind[(i,j)])
            objective_terms.append(-C*z[(i,j)]*data['widths'][i]*data['lengths'][i])
            #objective_terms.append(-D*y[(i,j)]*data['widths'][i]*data['heights'][i])
            #objective_terms.append(-E*x[(i,j)]*data['lengths'][i]*data['heights'][i])
    solver.Maximize(solver.Sum(objective_terms))

    status = solver.Solve()

    if status == pywraplp.Solver.OPTIMAL:
        print('Total objective value:', objective.Value())
        total_width = 0
        total_length = 0
        total_height = 0
        total_value=0
        cube_origins={}
        for j in data['bins']:
            bin_width = 0
            bin_length = 0
            bin_height = 0
            bin_value = 0
            cube_origins[j]={}
            cube_origins[j]['x']=[]
            cube_origins[j]['y']=[]
            cube_origins[j]['z']=[]
            print('Bin ', j, '\n')
            for i in data['items']:
                if ind[(i,j)].solution_value() > 0:
                    print('Item', i, '- width:', data['widths'][i], ' - length:', data['lengths'][i],' - height:',
                          data['heights'][i],' value:',data['values'][i])
                    bin_width += data['widths'][i]
                    bin_length += data['lengths'][i]
                    bin_height += data['heights'][i]
                    bin_value += data['values'][i]
                    print(x[(i,j)].solution_value())
                    print(y[(i,j)].solution_value())
                    print(z[(i,j)].solution_value())
                    cube_origins[j]['x'].append(x[(i,j)].solution_value())
                    cube_origins[j]['y'].append(y[(i,j)].solution_value())
                    cube_origins[j]['z'].append(z[(i,j)].solution_value())
                else:
                    cube_origins[j]['x'].append(-1)
                    cube_origins[j]['y'].append(-1)
                    cube_origins[j]['z'].append(-1)
                    print('Item', i, '- width:', data['widths'][i], ' - length:', data['lengths'][i],' - height:',
                          data['heights'][i],' value:',data['values'][i],'not in bin')

            print('Packed bin value:', bin_value)
            print()
            total_width += bin_width
            total_length += bin_length
            total_height += bin_height
            total_value += bin_value
            make_boxes(data,cube_origins,j)
        print('Total Packed Value: ', total_value)

    else:
        print(status)
        print('The problem does not have an optimal solution.')
# %%
if __name__ == '__main__':
    main()

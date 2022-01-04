# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 19:29:11 2021

@author: harres.tariq
"""

from ortools.linear_solver import pywraplp
#import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
# %%
def create_data_model():
    """Create the data for the example."""
    data = {}
    widths = [3,1,2,1,2,3]
    lengths = [2,2,3,1,2,1]
    values = [widths[i]*lengths[i] for i in range(len(widths))]#[1 for i in range(len(widths))]
    data['widths'] = widths
    data['lengths']=lengths
    data['values'] = values
    data['items'] = list(range(len(widths)))
    data['num_items'] = len(widths)
    num_bins = 2
    data['bins'] = list(range(num_bins))
    data['bin_max_width'] = [3,3]
    data['bin_max_length'] = [3,3]
    return data
# %%
def make_rectangles(data,cube_origins,binj):
    clr=['firebrick','yellow','green', 'red','blue','magenta','pink','cyan','khaki','gold','rosybrown','lightcoral','tan','olive']
    fig, ax = plt.subplots()
    ax.set_xlim([0,data['bin_max_width'][binj]])
    ax.set_ylim([0,data['bin_max_length'][binj]])
    for i in data['items']:
        try:
            if(cube_origins[binj]['x'][i]>=0):
                rect= patches.Rectangle((cube_origins[binj]['x'][i], cube_origins[binj]['y'][i]), data['widths'][i], data['lengths'][i], linewidth=5, facecolor=clr[i])
                ax.add_patch(rect)
        except:
            continue
        
    plt.show()
# %%
def main():
    data = create_data_model()

    # Create the mip solver with the SCIP backend.
    solver = pywraplp.Solver.CreateSolver('SCIP')

    x={}
    y={}    
    ind = {}
    for i in data['items']:
        for j in data['bins']:
            x[(i,j)]=solver.IntVar(0, data['bin_max_width'][j],'x'+str(i)+str(j))
            y[(i,j)]=solver.IntVar(0, data['bin_max_length'][j],'y'+str(i)+str(j))
            ind[(i, j)] = solver.IntVar(0, 1, 'ind_%i_%i' % (i, j))
    
    x_bool={}
    y_bool={}
    M=int(1e6)
    for i in data['items']:
        for j in data['bins']:
            for k in data['items']:
                if i>=k:
                    continue
                else:
                    x_bool[(i,j,k)] = solver.IntVar(0,1,"x_bool_"+str(i)+str(j)+str(k))
                    y_bool[(i,j,k)] = solver.IntVar(0,1,"y_bool_"+str(i)+str(j)+str(k))  
                    x_bool[(k,j,i)] = solver.IntVar(0,1,"x_bool_"+str(k)+str(j)+str(i))
                    y_bool[(k,j,i)] = solver.IntVar(0,1,"y_bool_"+str(k)+str(j)+str(i))
                    solver.Add(x[(i,j)]+data['widths'][i]<=x[(k,j)]+M*x_bool[(i,j,k)])
                    solver.Add(y[(i,j)]+data['lengths'][i]<=y[(k,j)]+M*y_bool[(i,j,k)])
                    solver.Add(x[(k,j)]+data['widths'][k]<=x[(i,j)]+M*x_bool[(k,j,i)])
                    solver.Add(y[(k,j)]+data['lengths'][k]<=y[(i,j)]+M*y_bool[(k,j,i)])
                    tmp=[]
                    tmp.append(x_bool[(i,j,k)])
                    tmp.append(y_bool[(i,j,k)])
                    tmp.append(x_bool[(k,j,i)])
                    tmp.append(y_bool[(k,j,i)])
                    solver.Add(sum(tmp)<=3)

        for j in data['bins']:
            solver.Add(x[(i,j)]+data['widths'][i]<=data['bin_max_width'][j]+M*(1-ind[(i,j)]))
            solver.Add(y[(i,j)]+data['lengths'][i]<=data['bin_max_length'][j]+M*(1-ind[(i,j)]))
    

    
    for i in data['items']:
        solver.Add(sum(ind[(i, j)] for j in data['bins']) <= 1)
    for j in data['bins']:    
        solver.Add(sum(ind[(i,j)]*data['values'][i] for i in data['items'])<=
               data['bin_max_width'][j]*data['bin_max_length'][j])
    

    # Objective
    objective = solver.Objective()

    for i in data['items']:
        for j in data['bins']:
            objective.SetCoefficient(ind[(i,j)], data['values'][i])
    objective.SetMaximization()

    status = solver.Solve()

    if status == pywraplp.Solver.OPTIMAL:
        print('Total packed value:', objective.Value())
        total_width = 0
        total_length = 0
        cube_origins={}
        for j in data['bins']:
            bin_width = 0
            bin_length=0
            bin_value = 0
            cube_origins[j]={}
            cube_origins[j]['x']=[]
            cube_origins[j]['y']=[]
            print('Bin ', j, '\n')
            for i in data['items']:
                if ind[(i,j)].solution_value() > 0:
                    print('Item', i, '- width:', data['widths'][i], ' - length:', data['lengths'][i],
                          ' value:',data['values'][i])
                    bin_width += data['widths'][i]
                    bin_length += data['lengths'][i]
                    bin_value += data['values'][i]
                    print(x[(i,j)].solution_value())
                    print(y[(i,j)].solution_value())
                    cube_origins[j]['x'].append(x[(i,j)].solution_value())
                    cube_origins[j]['y'].append(y[(i,j)].solution_value())
                    for k in data['items']:
                        if i>=k:
                            continue
                        else:
                            print('i,j,k: ',i,j,k)
                            print(x_bool[(i,j,k)].solution_value())
                            print(y_bool[(i,j,k)].solution_value())
                            print(x_bool[(k,j,i)].solution_value())
                            print(y_bool[(k,j,i)].solution_value())
                else:
                    cube_origins[j]['x'].append(-1)
                    cube_origins[j]['y'].append(-1)
                    print('Item', i, '- width:', data['widths'][i], ' - length:', data['lengths'][i],
                          ' value:',data['values'][i],'not in bin')

            print('Packed bin value:', bin_value)
            print()
            total_width += bin_width
            total_length += bin_length
            make_rectangles(data,cube_origins,j)
        print(cube_origins)
    else:
        print(status)
        print('The problem does not have an optimal solution.')
# %%
if __name__ == '__main__':
    main()

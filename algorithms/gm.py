import os
from math import *
import argparse
from collections import Counter
from numpy import *


def Gauss(dis,sigma2):
    factor1=1/sqrt(2*pi*sigma2)
    res=factor1*exp(-1*dis**2/2/sigma2)
    return res


def get_distance(lon1, lat1, lon2, lat2):
    # return distance (kilometer) of two GPS locations
    lon1 = float(lon1)
    lat1 = float(lat1)
    lon2 = float(lon2)
    lat2 = float(lat2)
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * asin(sqrt(a))
    r = 6378.137
    return c * r 


def gm_weight(loc_trace, isp_traces):
    #=================================Parameters of  GM===============================================
    Dmax=20
    p_i = [0.359328228655128,0.127758338037376,0.0447893488317361,0.0241707762264612,0.0297876415253099,0.0165717874252577,0.0216174941739561,0.0179831329718711,0.00954045567345170,0.0205600290567751,0.0221493117252212,0.0242575832039284,0.0325897392623172,0.0232985584468549,0.0159380782247211,0.0180647434662568,0.0340761131469641,0.0268338177931426,0.0111038976455218,0.0215709134021413,0.0151585156240235,0.0181553077403627,0.0378010615753326,0.0268951261658895,2.90759837461854e-181]
    sigma2 = [0.405744072706171,0.257764914697904,0.465629601281041,1.13993218413085,0.262573868741368,1.05455131004633,0.205009310987163,0.777525394649875,1.20241270632169,0.370100173232864,1.19870403406937,0.311848806748384,0.208668858134648,0.290757301745692,0.631637191135591,0.217271429723380,0.182220595450815,0.186990941668596,1.24754204640114,0.178123948703036,1.21203509798196,0.879854600932061,0.285474658908356,0.796374487384355,Dmax**2]
    Nh=25
    Ezero=1
    #==================================================================================================
    try:
        loc_points = dict(zip([k.strip().split(",")[0] for k in loc_trace.split(";")[1].split("|") if len(k.strip()) > 0],
                              [k.strip().split(",")[1] for k in loc_trace.split(";")[1].split("|") if len(k.strip()) > 0]))
    except IndexError:
        print 'loc_trace', loc_trace
        print 'isp_trace', isp_trace
        return float(-Inf)

    weight_list = []
    for isp_trace in isp_traces:
        isp_id = isp_trace.split(";")[0]
        try:
            isp_points = dict(zip([k.strip().split(",")[0] for k in isp_trace.split(";")[1].split("|")
                                   if len(k.strip()) > 0 and len(k.split(",")) >= 2],
                                  [k.strip().split(",")[1] for k in isp_trace.split(";")[1].split("|")
                                   if len(k.strip()) > 0 and len(k.split(",")) >= 2]))
                                   

            p_zero=p_i[Nh-1]*Gauss(Dmax,sigma2[Nh-1])

            Nlog=len(isp_points.values())
            LogHist=Counter(isp_points.values())
            Nrg=len(LogHist.keys())
            Edge=dict()
            for rg in LogHist.keys():
                Edge[rg]=float(Ezero+LogHist[rg])/(Nlog+Ezero*Nrg)
            
            
            weight=0
            for idx in loc_points.keys():
                p=p_zero
                locLOC=loc_points[idx]
                for h in range(0,Nh-1):
                    idh=str(int(idx)-h)
                    if isp_points.has_key(idh):
                        rg=isp_points[idh]
                        locRG = region_c_dict[rg]
                        dis = get_distance(float(locLOC.split("_")[1]), float(locLOC.split("_")[0]), float(locRG.split("_")[0]),
                                           float(locRG.split("_")[1]))
                        p=p+p_i[h]*Gauss(dis,sigma2[h])
                    else:
                        for rg in Edge.keys():
                            locRG = region_c_dict[rg]
                            dis = get_distance(float(locLOC.split("_")[1]), float(locLOC.split("_")[0]), float(locRG.split("_")[0]),
                                               float(locRG.split("_")[1]))
                            p=p+p_i[h]*Gauss(dis,sigma2[h])*Edge[rg]
                            
                weight=weight+log(p)
                                   
            weight_list.append((isp_id,weight))
        except IndexError:
            print 'data error'
    return weight_list


def hit_precision(rank, k):
    # return hit-precision with parameter K=k
    if rank > k:
        rank = k
    return 1 - rank * 1.0 / k


def single_task(locfile_name, input_path, output_path, save_weights=False):
    # deanonymization task of single user
    fr = open(input_path + locfile_name, 'r')
    all_traces = fr.readlines()
    loc_trace = all_traces[0]
    isp_traces = [all_traces[k] for k in range(len(all_traces)) if k >= 1]
    ispUID_score = gm_weight(loc_trace, isp_traces)
    ispUID_score_sorted = sorted(ispUID_score, key=lambda x: x[1], reverse=True)
    ispUID_sorted = [x[0] for x in ispUID_score_sorted]
    locUID = loc_trace.split(";")[0]
    if locUID in ispUID_sorted:
        rank_true = ispUID_sorted.index(locUID)
    else:
        rank_true = len(ispUID_sorted)
    print rank_true
    acc = hit_precision(rank_true, K_hp)
    if save_weights:
        fw = open(output_path + locfile_name, 'w')
        for i in range(len(ispUID_score_sorted)):
            fw.write(str(ispUID_score_sorted[i][0]) + '\t' + str(ispUID_score_sorted[i][1]) + '\n')
        fw.close()
    return locfile_name, acc


def main(input_path='./input/', output_path='./output/', save_weights=False):
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    file_names = os.listdir(input_path)
    acc_list = []
    for file_name in file_names:
        acc_list.append(single_task(file_name, input_path, output_path, save_weights))
    with open(output_path + '/result', 'w') as fw:
        for acc in acc_list:
            fw.write(str(acc[0]) + '\t' + str(acc[1]) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', nargs='?', default='./input/', help='Path of input data')
    parser.add_argument('--output_path', nargs='?', default='./output/gm/', help='Path of output data')
    parser.add_argument('--region_centre', nargs='?', default='./RegionCenters', help='Path of region centres\' location ')
    parser.add_argument('--save_weights', nargs='?', default=False, help='Save weight or not')
    parser.add_argument('--scale', nargs='?', default=0.01, help='Spatial scale', type=float)
    parser.add_argument('--K', nargs='?', default=10, help='K in hit-precision', type=int)
    parser.add_argument('--Tmax', nargs='?', default=24*8, help='The number of time-bins', type=int)

    args = parser.parse_args()
    K_hp = args.K
    region_centre = args.region_centre
    region_c_list = open(region_centre, 'r').readlines()
    region_c_dict = dict(
        zip([k.split("_")[0] for k in region_c_list], [k.split("_")[1] + "_" + k.split("_")[2] for k in region_c_list]))

    Tmax=args.Tmax
    AllTimeSlot=range(0,Tmax)

    main(input_path=args.input_path, output_path=args.output_path, save_weights=True)

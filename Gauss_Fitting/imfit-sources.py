from readcol import *
import numpy as np



def extractcomponents(sourcename=None,componentlist=None,comp=None,rms=None):
   print comp
  # outputfilefits = open(sourcename+'-fitresults.txt','w')
   #outputfileseps = open(sourcename+'-separations.txt','w')
   cl.open(componentlist)
   separations=np.zeros(comp, dtype=float)
   e_separations=np.zeros(comp, dtype=float)
   rakeep=np.zeros(comp, dtype=float)
   deckeep=np.zeros(comp, dtype=float)
   e_rakeep=np.zeros(comp, dtype=float)
   e_deckeep=np.zeros(comp, dtype=float)

   for j in range(comp):
      fit=cl.getcomponent(j)
      flux = fit['flux']['value']
      e_flux = fit['flux']['error']
      print flux[0]
      ra = fit['shape']['direction']['m0']['value']*180.0/3.14159265
      dec = fit['shape']['direction']['m1']['value']*180.0/3.14159265
      e_ra = fit['shape']['direction']['error']['latitude']['value']
      e_dec = fit['shape']['direction']['error']['longitude']['value']
      bmaj = fit['shape']['majoraxis']['value']*60.0     
      bmin = fit['shape']['minoraxis']['value']*60.0       
      pa =   fit['shape']['positionangle']['value']
      print ra,dec,e_ra,e_dec,bmaj,bmin
      line = str(sourcename[:])+" "+str(flux[0])+" "+str(e_flux[0])+ " " +str(ra)+" "+str(dec)+" "+str(e_ra)+" "+str(e_dec)+" " +str(bmaj)+" "+str(bmin)+" "+str(pa)+" "+str(rms)+"\n"
      rakeep[j]=ra
      deckeep[j]=dec
      e_rakeep[j]=e_ra
      e_deckeep[j]=e_dec
      outputfilefits.writelines(line)
   print range(comp)
   iter=0
   for j in range(comp):
      for k in range(comp):
         print e_rakeep[j], e_deckeep[j],e_rakeep[k], e_deckeep[k]
         if j != k and k > j:
            print iter
            dra=rakeep[j]* 3600.0-rakeep[k]* 3600.0
            ddec=deckeep[j]* 3600.0-deckeep[k]* 3600.0

            e_dra=(e_rakeep[j]**2 + e_rakeep[k]**2)**0.5
            e_ddec=(e_deckeep[j]**2 + e_deckeep[k]**2)**0.5

            separations[iter]=sqrt( (dra)**2 + (ddec)**2)  
            e_rasq=(dra*e_dra)**2
            e_decsq=(ddec*e_ddec)**2

            print e_rasq, e_decsq
            e_separations[iter] =  sqrt( e_rasq + e_decsq )/separations[iter]
            iter += 1

   #for j in range(iter):
     #line = str(separations[j])+ " "+str(e_separations[j])+ " "+str(separations[j]*dist) + " "+str(e_separations[j]*dist) +"\n"
     # outputfileseps.writelines(line)
   cl.close()
#   outputfilefits.close()
#   outputfileseps.close()
   print 'end of extract'


dist=230.0

#sources=['Per2','Per11','Per12','Per17','Per18','Per22','Per26','Per27','Per33','Per35','Per36','Per40','Per44','Per48','Per49','Per55','Per106','Per107','Per120','Per122']

#comps=[2,3,2,2,3,2,2,2,4,2,2,2,5,2,2,2,2,2,2,2]
inp(taskname='imfit')
cl.close()
imagelisttemp=readcol("images.txt")
#imagelist=imagelisttemp.array().tolist()
imagelist=[]
for image in imagelisttemp:
   imagelist.append(image[:][0])
#   print image

#print imagelist

sourcestemp=[line.replace('.B.Ka.cont.image.tt0','') for line in imagelist]
sources=[source.replace('../','') for source in sourcestemp]
#print imagelist
#print sources

sources.append('Per110')
sources.append('Per21')
sources.append('SVS13B')
sources.append('IRAS4C')

imagelist.append('../Per110.B.Ka.cont.image.tt0')
imagelist.append('../Per21.B.Ka.cont.image.tt0')
imagelist.append('../SVS13B.B.Ka.cont.image.tt0')
imagelist.append('../IRAS4C.B.Ka.cont.image.tt0')

outputfilefits = open('all-fitresults.txt','w')

#sources=['Per1','Per2']
print len(sources)

for i in range(len(sources)):
  print i
  if sources[i] != 'Per101' and sources[i] != 'Per103' and sources[i] != 'Per117' and sources[i] != 'Per4' and sources[i] != 'Per39' and sources[i] != 'Per43' and sources[i] != 'Per45' and sources[i] != 'Per59' and sources[i] != 'Per60':
   inp(taskname='imfit')
   region='../regions-B/'+sources[i]+'-B.rgn'
   estimates=sources[i]+'-gaussest.txt'
   imagename=imagelist[i]
   complist=sources[i]+'-fitcomps.cl'
   print sources[i]
   overwrite=True   
   inp(taskname='imfit')
   imfit()
   ia.open(imagename)
   stats=ia.statistics()
   rms=stats["rms"][0]
   ia.close()
   extractcomponents(sources[i],complist,1,rms)
   
      


outputfilefits.close()




ap.experiment <- function(data , r=6 , window.size=1000 , slide.size=10, q=0.01)
{
  protmime <- 0
  m3 <- 0
  h=0
  xx <- matrix(0,1,210)
  xxx <-1
 # m1 <-memory.size(T)
 #gc1 <- gc()
 # print(gc1)
  for(i in seq(1,(209*slide.size+1),slide.size))
    
    #for(i in seq(1,(29*slide.size+1),slide.size))
  #for(i in c(10,110,210,310,400,510,620,730,800,900))
  {
  #i=1
  #for 循环
  data.1 <- data[1:(1+window.size-1),]
 data.2 <- data[c(1:(1+window.size-1),(10000+i):(10000+i-1+slide.size)),]
 # data.1 <- data[i:(i+window.size-1),]
 #data.2 <- data[i:(i+window.size-1+slide.size),]
  apres.1 <- apcluster(negDistMat(r=r),data.1,q=q)
  apres.2 <- apcluster(negDistMat(r=r),data.2,q=q)
  exemplars.1 <- as.vector(apres.1@exemplars) 
  exemplars.2 <- as.vector(apres.2@exemplars)
  #exemplars.1 <- apres.1@exemplars
  #exemplars.2 <- apres.2@exemplars
  print('一共聚了多少类')
  print(length(exemplars.2))
  tar <- 0
  #timestart<-Sys.time()
  DDD1<- 1000
  #print(apres.1)
 #print(apres.2)
#  print(memory.size(max = F))
  #print(exemplars.1)
 # print(exemplars.2)
 # ptm <- proc.time()
  
  rm(con.table)
  rm(c1)
  rm(E.nij)
  rm(x2)
  rm(ap.dist)
  rm(SOE)
  rm(D)
  rm(D1)
  rm(D2)
  rm(con)
  rm(SOE.all)
  rm(weight)
  gc(reset = T)
  
  
  
  
  
############################################################  
  
 ap.dist <- as.matrix(dist(data.2))
 ap.dist1  <<-  ap.dist
 D <- vector(mode = 'integer')
 con <-vector()
 SOE <-vector()
for(j in 1:length(exemplars.2))
  {
    D <- as.vector(unlist(apres.2@clusters[j]))
    D1 <- D[D<=(window.size)]
   D2 <- D[D>(window.size)]
    
    D1 <- D[D<=(window.size)]
    
     D2 <- D[D>(window.size)]
      
  if (length(D2)!=0&&length(D1)!=0)
    {
     D1 <- mean(ap.dist[exemplars.2[j],D1])
      D2 <- mean(ap.dist[exemplars.2[j],D2])
      if(D1<DDD1)
        DDD1 <- D1
      con <- D1/D2
     
    }
    
  }
  
  
  ###############################################
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
 # print(memory.size())
  #m1 <- memory.size()
  if (length(exemplars.1)!=length(exemplars.2))
     { print("聚类数目发生了改变")
      tar <- 1 
      xx[1,xxx]<- 10
     }
     else
     {
    #  print(exemplars.2)
     # print(exemplars.2[length(exemplars.2)])
     #  print(i+window.size-1)
       if(exemplars.2[length(exemplars.2)]+i > (i+window.size-1))
       {
         print("exemplar在active中，发现离群点")
         tar <- 1
         xx[1,xxx]<- 11
       }
         else
         {
         #  print('11111111')
           #计算CT
           con.table <- matrix(nrow = length(exemplars.1)+1,ncol = length(exemplars.1)+1)
         # print('111111111')
            #  print(con.table)
           #  print(exemplars.1)
           for(i in 1:length(exemplars.1))
           {
             for(j in 1:length(exemplars.1))
             {
               #intersect可以比较出2个向量中相同部分的个数
               con.table[i,j] <-length(intersect(as.vector(unlist(apres.1@clusters[i])),as.vector(unlist(apres.2@clusters[j]))))
               
             }
           }
        #   print(con.table)
            c1 <- max.col(con.table[1:length(exemplars.1),1:length(exemplars.1)])
           
           for(i in 1:length(exemplars.1))    # set max to 0
           {
             
             #con.table[i,(max.col(con.table)[i])] <-0
             #con.table[i,max.col(con.table)] <-0
             con.table[i,c1[i]] <- 0
           }
           #con.table.test <<- con.table
         #  print(con.table)
           for(i in 1:length(exemplars.1))  # compute nii
           {
             
             con.table[i,(length(exemplars.1)+1)]<- sum(con.table[i,1:length(exemplars.1)])
             con.table[(length(exemplars.1)+1),i]<- sum(con.table[1:length(exemplars.1),i])
             
           }
           
        #  print("1111111111")
        #   print(con.table)
           con.table[(length(exemplars.1)+1),(length(exemplars.1)+1)] <- sum(con.table[1:length(exemplars.1),(length(exemplars.1)+1)])  #compute n ,n =con.table[(exemplars.1+1),(exemplars.1+1)]
           
           
         # print("1111111111")
          # print(con.table)
           
           #identical可以比较两个向量是否相同
      # if (!identical(exemplars.1,exemplars.2))  #如果两个聚类中心点不相同
         
     ####  {  
         
       
        x2 <- 0
       # print(x2)
        if (con.table[(length(exemplars.1)+1),(length(exemplars.1)+1)]!=0)
        {
       for(i in 1:length(exemplars.1))   #compute x2
       {
         for(j in 1:length(exemplars.1))
         {
           
           E.nij <- con.table[i,(length(exemplars.1)+1)]*con.table[(length(exemplars.1)+1),j]/con.table[(length(exemplars.1)+1),(length(exemplars.1)+1)]
          if (E.nij!=0)
            x2 <-( x2 + ((con.table[i,j] - E.nij)^2 / E.nij ))
        #  print("EIJ:") 
        #  print(E.nij)
        #  print('X2:')
        #   print(x2)
         }
       }
        }
     #   print('x2总')
      #  print(x2)
       if (x2 > 3.84)
       {
         print('95%发现异常')
         tar <- 1
         xx[1,xxx]<- 12
       }
        
    #####   #}
           #else  #两个聚类中心点相同。
        if (tar !=1)
           {
             
             
             
             ap.dist <- as.matrix(dist(data.2))
            # ap.dist1  <<-  ap.dist
            # D <- vector(mode = 'integer')
             con <-vector()
             SOE <-vector()
             for(j in 1:length(exemplars.2))
             {
               D <- as.vector(unlist(apres.2@clusters[j]))
               D1 <- D[D<=(window.size)]
              D2 <- D[D>(window.size)]
            
               #  D1 <- D[D<=(window.size)]
            
               #   D2 <- D[D>(window.size)]
        #      print('D1:')
       #     print(D1)
       #     print('D2:')
      #       print(D2)
             # print('D:')
              #print(D)
               if (length(D2)!=0&&length(D1)!=0)
               {
               D1 <- mean(ap.dist[exemplars.2[j],D1])
               D2 <- mean(ap.dist[exemplars.2[j],D2])
               if(D1<DDD1)
               DDD1 <- D1
               con <- D1/D2
               }
               else 
                 con <- 0
              
              
              if(con.table[(length(exemplars.2)+1),(length(exemplars.2)+1)]!=0)
                weight<-con.table[(length(exemplars.2)+1),j]/con.table[(length(exemplars.2)+1),(length(exemplars.2)+1)]
               else
                 weight <- 1/length(exemplars.2)
                SOE[j] <- 2*con/(1+con*weight)
            # print('con是')
           #  print(con)
           #  print('weight是')
           #  print(weight)
                
                }
            #  print('SOE是')
            # print(SOE)
              SOE.all <- min(SOE[SOE>0])
              
               print("SOE的值为：") 
               print(SOE.all)
               xx[1,xxx]<- SOE.all
               
               print('阈值是：')
               print((2*(1/2))/(1+(1/2)*(1/length(exemplars.2))))
          #    print((2*(2/3))/(1+(2/3)*(1/length(exemplars.2))))
             #  print(length(exemplars.2))
             #  if (SOE.all < (2*0.5)/(1+0.5*(1/length(exemplars.2))))
                 if (SOE.all < (2*(1/2))/(1+(1/2)*(1/length(exemplars.2))))
               #    if (SOE.all < (2*(2/3))/(1+(2/3)*(1/length(exemplars.2))))
              {
                print("SOE检测出异常")
                
                tar <- 1
               }
             
           }
           
           
           
           
           
           
         }
     }

    #else
  #{
 #  两个聚类中心点编号相同 
 #  }
  if (tar==0)
    print('未发现异常')
 # print('DDD1:')
 # print(DDD1)
 print("次数：")
  h=h+1
  print(h)
  print('     ')
  
  
 #  gc3 <- gc()
  # print(gc3)
  #print(apres.1)
  #print(apres.2)
  
  
  #return(con.table)
  
#  timeend<-Sys.time()
 # runningtime<-timeend-timestart
  #print(runningtime)
  
  
  #proc.time() - ptm
 # print(proc.time() - ptm)
 # protmime <- protmime+(proc.time() - ptm)
 # print('内存是')
  
 # print(memory.size(F))
 # m2 <-memory.size(F)
 # m3 <- m2-m1
 # print(m3)
# print(memory.size(max = F))
  #m3 <-m3+m2-m1
 # print(m2-m1)
  
 # print(ptm)
  
  xxx <- xxx+1
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}
  
  #slide10.window2000.vv2.9.30 <<-xx
  slide10.window2000.vv2.4.12 <<-xx
#  m2 <- memory.size(T)
 # gc2 <- gc()
 # print(gc2)
  #print(protmime/10)
#  print('内存!!!!是')
 # total <- (as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", system("wmic OS get FreePhysicalMemory /Value", intern=TRUE))[3]))/1000
 # print(total)
 # print(m3/10)
 # print(m3/10)
# print(m2-m1)
# print(gc2[4]-gc1[4])
  
}
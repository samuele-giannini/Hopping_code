      program morphology
      real T(3),U(3,3),C(4,3) 
      real R(8000,3)
      integer L(10,6)

      open (unit=30,file='coordinates.dat')
      open (unit=31,file='connectivity.dat')
      open (unit=32,file='pbc.dat')
       
c generate morphology.prm (fort.31), coord.xyz (fort.30), 
c pbc.inp (fort.33) to be used by the SCD suite of programs
      open (unit=10,file='gm.inp')
      read (10,*)nmuc
c read number of molecules in the unit cell (nmuc) and their
c coordinates C(k,j)
      do k=1,nmuc
       read (10,*)(C(k,j),j=1,3)
      enddo
      read (10,*)
c read unit cell vector U(m,j)
      do m=1,3
       read (10,*)(U(m,j),j=1,3)
      enddo
      read (10,*)
c read supecell size to be built
      read (10,*)na,nb,nc
      read (10,*)
c read number of unique interactions to be considered
      read (10,*)ni
      read (10,*)
      do i=1,ni
       read (10,*)(L(i,k),k=1,6)
       write(*,*)(L(i,k),k=1,6)
      enddo
      write (31,*)na*nb*nc*nmuc
c======= cycle over the supercell
      do ia=0,na-1
      do ib=0,nb-1
      do ic=0,nc-1
c       translational vector for this supercell
        T(1)=ia*U(1,1)+ib*U(2,1)+ic*U(3,1)
        T(2)=ia*U(1,2)+ib*U(2,2)+ic*U(3,2)
        T(3)=ia*U(1,3)+ib*U(2,3)+ic*U(3,3)
c       write the coordinates
        do i=1,nmuc

c        is1 should grow continuously
         is1=nmuc*(nb*nc*ia+nc*ib+ic)+i
         do k=1,3
          R(is1,k)=C(i,k)+T(k)
         enddo
         write (30,30)(R(is1,k),k=1,3)

        enddo
c the molecule m in the replica ia,ib,ic is molecule 
c number (na*nb*ia+nb*ib+ic)+m

c       write the interactions
        do i=1,ni
c the first molecule of the interacting pair is
         is1=nmuc*(nb*nc*ia+nc*ib+ic)+L(i,1)
c to find which one is the second molecule
         ia2=ia+L(i,3)
         if (ia2.eq.-1) then 
          ia2=na-1
         endif
         if (ia2.eq.na) then 
          ia2=0
         endif
         ib2=ib+L(i,4)
         if (ib2.eq.-1) then 
           ib2=nb-1
         endif
         if (ib2.eq.nb) then 
           ib2=0
         endif
         ic2=ic+L(i,5)
         if (ic2.eq.-1) then 
           ic2=nc-1
         endif
         if (ic2.eq.nc) then 
           ic2=0
         endif
         is2=nmuc*(nb*nc*ia2+nc*ib2+ic2)+L(i,2)
c write the interacting pair
         write (31,31)is1,is2,L(i,6)
        enddo
      enddo
      enddo
      enddo
c=======
c     fort.33 will be the pbc.inp file
      write (32,30)(na*U(1,j),j=1,3)
      write (32,30)(nb*U(2,j),j=1,3)
      write (32,30)(nc*U(3,j),j=1,3)
     
      stop
  30  format (3F10.4,I6)
  31  format (2I6,I7,3I6)
  32  format (3F10.4,I8) 
      stop
      end


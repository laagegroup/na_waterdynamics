PROGRAM integtimes

USE calc

implicit none

integer              ::s,n,i,ib,ntcf,l,nout,nfit,nsum,t
parameter(ntcf=19,nout=20,nfit =21,nsum=22)
!input
integer              ::site_number,stat,nlines,integmax
logical              ::tcf_exist
character*4          ::site
character*2          ::blk
real*8,dimension(:),allocatable   ::time
real*8,dimension(:),allocatable   ::tcf
logical   ::decay,noise,iffit
real*8                            ::dt
real*8                            ::integt,moyt,moyt2
real*8,dimension(:),allocatable   ::integtime
real*8                            ::stdev
real*4                            ::noise_max
integer                           :: nblock
!fit
integer                                    ::length_fit
integer                                    ::dimfit,itfitmin,itfitmax
real*8                                     ::a,b,siga,sigb,chi2,q
real*8,dimension(:),allocatable            ::xfit,yfit,sig
!correction
real*8                                     ::corrint

print *,"number of sites ?"
read(*,*) site_number
print *,"number of blocks?"
read(*,*) nblock
print *,"fit length? (ps)"
read(*,*) length_fit  !fit over the last xxx ps

noise_max=0.2  
!print *,"noise max on lnC(t) ?"
!read(*,*) noise_max  !fit over the last xxx ps


print *,'nsite =',site_number
allocate(integtime(site_number+1))
integtime=0
print *,'lengthfit =',length_fit
print *,'maxnoise =', noise_max

nlines = 0

open(nsum,file = 'summary_integtime_jump.dat')

do s = 1,site_number !for each site
    tcf_exist = .false.
    site = intchar4(s)
    inquire(file='jcf.out'//site, exist=tcf_exist) !tests if file exists 
    print *,s

    if (tcf_exist) then

        !analysis of total tcf
        nlines=0
        if (nlines.eq. 0) then !first time : counts the number of lines in file to allocate time and tcf arrays
            stat=0
            open(ntcf,file = 'jcf.out'//site, readonly, status = 'old')
            do while(stat >=0)
                read(ntcf,*,iostat = stat)
                nlines = nlines +1
            end do
            close(ntcf)
            allocate(time(nlines),tcf(nlines))
        end if

        open(ntcf,file = 'jcf.out'//site,readonly,status='old')
           stat = 0
           l=0
           do while (stat >=0) !while file isnot finished
                l = l+1
                read(ntcf,fmt='(2f20.10)',iostat=stat) time(l), tcf(l) 
           end do
        close(ntcf)
        print *, 'end reading file jcf.out'//site
    
        !calculate integral of C(t) and correction if necessary
        ! integrate until goes to 0 or until the end or noise too high 
        decay=.false.
        noise=.false.
        do i=1,nlines-2
            ! if not end of decay
            if ((tcf(i+1).gt.0).and.(.not.(decay)).and.(.not.(noise))) then
                itfitmax= i
                integtime(s) = integtime(s) + (tcf(i+1) + tcf(i))/2d0*(time(i+1) - time(i)) 
            end if
            if ((tcf(i+1).le.0)) then
                decay=.true.
            end if
            if ((i.gt.10).and.( abs(log(tcf(i+1)) - log(tcf(i))) .gt. noise_max)) then
                noise=.true.
            end if
        end do

        dt=time(2) - time(1)

        open(nout,file='integ_jcf.out'//site)
        write(nout,*) 'integral until ', itfitmax*dt, '=', integtime(s)
        write(nout,*) 'maxnoise =  ', noise_max

        !linear fit on the last length_fit (ps) before end of tcf    
        !if not fully decayed:
        if ((.not.decay).or.(noise)) then
                dimfit = int(float(length_fit) / dt) +1 !number of points to fit
                itfitmin = itfitmax - dimfit +1 
                !need to fit on at least 5 ps
                iffit=.true.
                if (itfitmax*dt .le.8) then
                        write(nout,*) 'itfitmax too small'
                        print *, 'itfitmax too small'
                        iffit=.false.
                else if (itfitmin*dt .le.3) then
                        print *,'min too small'
                        itfitmin = int(3.0/dt)
                        dimfit = itfitmax - itfitmin +1
                end if
        !print *,"itfitmin=",itfitmin,"itfitmax=",itfitmax,'dimfit =',dimfit
        
                if (iffit.eqv..true.) then
                        write(nout,*) 'fit in the SHELL y = a + bx'
                        write(nout,*) 'itfitmin = ',itfitmin*dt, 'itfitmax = ',itfitmax*dt

                        allocate(xfit(dimfit),yfit(dimfit),sig(dimfit))
                        do i=1,dimfit
                                xfit(i) = time(itfitmin +i -1)
                                yfit(i) = log(tcf(itfitmin +i -1))
                                sig(i) = 1.0
                        end do

                        call fit(xfit,yfit,dimfit,sig,1,a,b,siga,sigb,chi2,q)
   
                        write(nout,*) 'with a=',a,'and b=',b
                        write(nout,*) 'siga=',siga,'and sigb=',sigb
                        write(nout,"('chi2 =',f8.3)") chi2
                        write(nout,"('q =',f8.3)") q
                        if (b.le.(-1d-10)) then
                                write(nout,"('taulong =',f12.3,' +/- ',f12.3)") -1d0/b,sigb/(b**2)
                        end if

                        deallocate(xfit,yfit,sig)
    
                        !calculate correction to the integral based on the fit
                        if (b.le.(-1d-10)) then
                                ! print *,'be lower than 10-6'
                                corrint=dexp(a)*(-1/b)*dexp(b*real(itfitmax)*dt)
                                integtime(s) = integtime(s) + corrint
                                write(nout,"('correction to the integral =',f12.4)") corrint
                        else 
                                write(nout,*) 'ERROR b>0 !!!'
                        end if
        

                        write(nout,"('integrated orcorr time after correction in the shell=',f12.3)") integtime(s)
                        !close(nout) 

                        !write the fit function in separated file
                        !open(nfit,file='fit.out'//site,form='formatted')
                        !write(nout,*) "# time(ps), ln(jcf), fit"
                        !do t=1,itfitmax
                        !        write(nfit,fmt='(f15.3,f15.7,f15.7)') real(t-1)*dt, log(tcf(t)), a+b*dble(t-1)*dt
                        !end do
                        !close(nfit)
                end if !iffit
        end if
        deallocate(time,tcf)


    ! repeat integration for all the blocks and calculate std deviation on
    ! integrated time over the n blocks.
        write(nout,*) "integtime block"
        moyt=0.0
        moyt2=0.0
        n=0.0
        stdev=0.0
        do ib=1,nblock
           blk = intchar2(ib)
           print *,site, "    ", blk
           inquire(file='jcf_block'//blk//'site'//site//'.out',exist=tcf_exist) 

           if (tcf_exist) then
                open(ntcf,file='jcf_block'//blk//'site'//site//'.out',readonly,status='old') 

                !calculation of integrated time on each block
                stat=0
                do while(stat >=0)
                        read(ntcf,*,iostat = stat)
                        nlines = nlines +1
                end do
                close(ntcf)
                allocate(time(nlines),tcf(nlines))

                open(ntcf,file='jcf_block'//blk//'site'//site//'.out',readonly,status='old') 
                stat = 0
                l=0
                do while (stat >=0) !while file isnot finished
                        l = l+1
                        read(ntcf,fmt='(2f20.10)',iostat=stat) time(l), tcf(l) 
                end do
                close(ntcf)
    
                !calculate integral of C(t) and correction if necessary
               ! integrate until goes to 0 or until the end or noise too high 
                decay=.false.
                noise=.false.
                integt=0.0
                corrint=0.0
                do i=1,nlines-2
                        ! if not end of decay
                        if ((tcf(i+1).gt.0).and.(.not.(decay)).and.(.not.(noise))) then
                                itfitmax= i
                                integt = integt + (tcf(i+1) + tcf(i))/2d0*(time(i+1) - time(i)) 
                        end if
                        if ((tcf(i+1).le.0)) then
                                decay=.true.
                        end if
                        if ((i.gt.10).and.( abs(log(tcf(i+1)) - log(tcf(i))) .gt. noise_max)) then
                                noise=.true.
                        end if
                end do

                !linear fit on the last length_fit (ps) before end of tcf    
                !if not fully decayed:
                if ((.not.decay).or.(noise)) then
                        iffit=.true.
                        dimfit = int(float(length_fit) / dt) +1 !number of points to fit
                        itfitmin = itfitmax - dimfit +1 
                        !need to fit on at least 5 ps
                        if (itfitmax*dt .le.8) then
                                write(nout,*) 'itfitmax too small'
                                print *, 'itfitmax too small'
                                iffit=.false.
                        end if
                        if (itfitmin*dt .le.3) then
                                print *,'min too small'
                                itfitmin = int(3.0/dt)
                                dimfit = itfitmax - itfitmin +1
                        end if
        
                        if (iffit.eqv..true.) then
                                allocate(xfit(dimfit),yfit(dimfit),sig(dimfit))
                                do i=1,dimfit
                                        xfit(i) = time(itfitmin +i -1)
                                        yfit(i) = log(tcf(itfitmin +i -1))
                                        sig(i) = 1.0
                                end do

                                call fit(xfit,yfit,dimfit,sig,1,a,b,siga,sigb,chi2,q)
                                deallocate(xfit,yfit,sig)
    
                                !calculate correction to the integral based on the fit
                                if ((b.le.(-1d-10)).and.(a.le.10)) then
                                        corrint=dexp(a)*(-1/b)*dexp(b*real(itfitmax)*dt)
                                        integt = integt + corrint
                                end if
                        end if !iffit
        
                end if

                deallocate(time,tcf)

                write(nout,*) integt
                moyt=moyt+integt
                moyt2=moyt2+integt**2
                n=n+1.0
             end if !tcf_exist
        end do !blocks

        print *,"moyt=",moyt,"moyt2=",moyt2,"n=",n
        if (n.gt.0) then
                moyt=moyt/float(n)
                moyt2=moyt2/float(n)
                print *,moyt2-moyt**2
                stdev=dsqrt(moyt2-moyt**2)
        end if

        write(nsum,*) s, integtime(s),stdev !en ps

    end if
    close(nout) 
    
end do !sites

close(nsum)

END PROGRAM integtimes

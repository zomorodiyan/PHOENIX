!______________________________________________________________________________
!
module entotemp
!______________________________________________________________________________
!
	use initialization
	use parameters
	use dimensions
	use species, only: concentration, mix, &
		tsolid2, tliquid2, acpa2, acpb2, acpl2

	contains
subroutine enthalpy_to_temp(ilo, ihi, jlo, jhi, klo, khi)
	implicit none
	integer, intent(in) :: ilo, ihi, jlo, jhi, klo, khi
	integer i,j,k
	! Composition-dependent locals
	real(wp) C_local, tsolid_l, tliquid_l, acpa_l, acpb_l, acpl_l
	real(wp) hsmelt_l, hlcal_l, deltemp_l, cpavg_l

	do k=max(1,klo-1),min(nk,khi+1)
!$OMP PARALLEL PRIVATE(C_local, tsolid_l, tliquid_l, acpa_l, acpb_l, acpl_l, &
!$OMP   hsmelt_l, hlcal_l, deltemp_l, cpavg_l)
!$OMP DO
	do j=max(1,jlo-1),min(nj,jhi+1)
	do i=max(1,ilo-1),min(ni,ihi+1)
		if (species_flag == 1) then
			C_local   = concentration(i,j,k)
			tsolid_l  = mix(tsolid,  tsolid2,  C_local)
			tliquid_l = mix(tliquid, tliquid2, C_local)
			acpa_l    = mix(acpa,    acpa2,    C_local)
			acpb_l    = mix(acpb,    acpb2,    C_local)
			acpl_l    = mix(acpl,    acpl2,    C_local)
			hsmelt_l  = acpa_l*tsolid_l**2/2.0_wp + acpb_l*tsolid_l
			cpavg_l   = (acpa_l*tsolid_l + acpb_l + acpl_l)*0.5_wp
			deltemp_l = max(tliquid_l - tsolid_l, 1.0_wp)
			hlcal_l   = hsmelt_l + cpavg_l*deltemp_l
		else
			tsolid_l  = tsolid;  tliquid_l = tliquid
			acpa_l    = acpa;    acpb_l    = acpb;    acpl_l = acpl
			hsmelt_l  = hsmelt;  hlcal_l   = hlcal;   deltemp_l = deltemp
		endif

		if(enthalpy(i,j,k).ge.hlcal_l) then
			fracl(i,j,k)=1.0
			temp(i,j,k)=(enthalpy(i,j,k)-hlcal_l)/acpl_l+tliquid_l
		elseif(enthalpy(i,j,k).le.hsmelt_l) then
			fracl(i,j,k)=0.0
			temp(i,j,k)=(sqrt(acpb_l**2+2*acpa_l*enthalpy(i,j,k))-acpb_l)/acpa_l
		else
			fracl(i,j,k)=(enthalpy(i,j,k)-hsmelt_l)/(hlcal_l-hsmelt_l)
			temp(i,j,k)=deltemp_l*fracl(i,j,k)+tsolid_l
		endif
	enddo
	enddo
!$OMP END PARALLEL
	enddo
	return

end subroutine enthalpy_to_temp

end module entotemp

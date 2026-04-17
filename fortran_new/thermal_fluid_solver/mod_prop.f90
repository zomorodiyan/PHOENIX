!______________________________________________________________________________
!
module property
!______________________________________________________________________________
!
	use initialization
	use parameters
	use dimensions
	use species, only: concentration, mix, &
		dens2, denl2, viscos2, tsolid2, tliquid2, &
		acpa2, acpb2, acpl2, thconsa2, thconsb2, thconl2, &
		pden2, pcpa2, pcpb2, pthcona2, pthconb2

	implicit none

	contains

!********************************************************************
subroutine properties(ilo, ihi, jlo, jhi, klo, khi)
	implicit none
	integer, intent(in) :: ilo, ihi, jlo, jhi, klo, khi
	integer i,j,k, jind
	real(wp) diffs, diffl
	real(wp) dL_mix, vMag, visT, diffT
	! Composition-dependent locals
	real(wp) C_local, tsolid_l, tliquid_l
	real(wp) dens_l, denl_l, viscos_l, acpa_l, acpb_l, acpl_l
	real(wp) thconsa_l, thconsb_l, thconl_l
	real(wp) pden_l, pcpa_l, pcpb_l, pthcona_l, pthconb_l

	do k=klo,khi
!$OMP PARALLEL PRIVATE(jind, dL_mix, vMag, visT, diffT, diffs, diffl, &
!$OMP   C_local, tsolid_l, tliquid_l, dens_l, denl_l, viscos_l, &
!$OMP   acpa_l, acpb_l, acpl_l, thconsa_l, thconsb_l, thconl_l, &
!$OMP   pden_l, pcpa_l, pcpb_l, pthcona_l, pthconb_l)
!$OMP DO
	do j=jlo,jhi
	do i=ilo,ihi

		if (species_flag == 1) then
			C_local   = concentration(i,j,k)
			tsolid_l  = mix(tsolid,  tsolid2,  C_local)
			tliquid_l = mix(tliquid, tliquid2, C_local)
			dens_l    = mix(dens,    dens2,    C_local)
			denl_l    = mix(denl,    denl2,    C_local)
			viscos_l  = mix(viscos,  viscos2,  C_local)
			acpa_l    = mix(acpa,    acpa2,    C_local)
			acpb_l    = mix(acpb,    acpb2,    C_local)
			acpl_l    = mix(acpl,    acpl2,    C_local)
			thconsa_l = mix(thconsa, thconsa2, C_local)
			thconsb_l = mix(thconsb, thconsb2, C_local)
			thconl_l  = mix(thconl,  thconl2,  C_local)
			pden_l    = mix(pden,    pden2,    C_local)
			pcpa_l    = mix(pcpa,    pcpa2,    C_local)
			pcpb_l    = mix(pcpb,    pcpb2,    C_local)
			pthcona_l = mix(pthcona, pthcona2, C_local)
			pthconb_l = mix(pthconb, pthconb2, C_local)
		else
			tsolid_l  = tsolid;  tliquid_l = tliquid
			dens_l    = dens;    denl_l    = denl
			viscos_l  = viscos;  acpa_l    = acpa
			acpb_l    = acpb;    acpl_l    = acpl
			thconsa_l = thconsa; thconsb_l = thconsb
			thconl_l  = thconl
			pden_l    = pden;    pcpa_l    = pcpa
			pcpb_l    = pcpb;    pthcona_l = pthcona
			pthconb_l = pthconb
		endif

		visT=0
		diffT=visT/0.9

		diffs=(thconsa_l*temp(i,j,k)+thconsb_l)/(acpa_l*temp(i,j,k)+acpb_l)
		diffl=thconl_l/acpl_l
		vis(i,j,k)=(viscos_l+visT)
		diff(i,j,k)=(diffl+diffT)
		den(i,j,k)=denl_l

		if(temp(i,j,k).ge.tliquid_l) cycle
			diff(i,j,k)=diffs
			vis(i,j,k)=vis_solid
			den(i,j,k)=dens_l

			if(z(nk)-z(k) .le. layerheight .and. solidfield(i,j,k) .le. powder_threshold)then
				den(i,j,k)=pden_l
				vis(i,j,k)=vis_solid
				diff(i,j,k)=(pthcona_l*temp(i,j,k)+pthconb_l)/(pcpa_l*temp(i,j,k)+pcpb_l)
			endif

		if(temp(i,j,k).le.tsolid_l) cycle
			diff(i,j,k)=(fracl(i,j,k)*diffl+(1.0-fracl(i,j,k))*diffs)
			vis(i,j,k)=(viscos_l+visT)
			den(i,j,k)=(fracl(i,j,k)*denl_l+(1.0-fracl(i,j,k))*dens_l)

	enddo
	enddo
!$OMP END PARALLEL
	enddo
end subroutine properties
end module property

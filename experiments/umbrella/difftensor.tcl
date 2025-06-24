cv colvar phi set collect_gradient 1
cv colvar psi set collect_gradient 1


# inverse molar masses (mol*g^-1)
variable invM {0.0832639  0.0713776  0.0832639  0.0832639  0.0713776}
# corresponding ids (4,6,8,14,16)


proc calc_difftensor { args } {
    cv colvar phi update
    cv colvar psi update
    set gradphi [cv colvar phi getgradients]
     # (4,6,8,14)
    set gradpsi [cv colvar psi getgradients]
    # (6,8,14,16)

    set dphi2 0.0
    set dpsi2 0.0
    set dphidpsi 0.0

    set phiix 0
    set psiix 0

    # compute contributions of shared atoms (ids 6,8,14)
    for {set i 1} {$i < 4} {incr i} {
        #puts [lindex $gradphi $i]
        #puts [lindex $gradpsi [expr {$i-1}]]
        foreach a [lindex $gradphi $i] b [lindex $gradpsi [expr {$i-1}]] {
                    set dphidpsi [expr {$dphidpsi + $a*$b*[lindex $::invM $i]}]
                    set dphi2 [expr {$dphi2 + $a*$a*[lindex $::invM $i]}]
                    set dpsi2 [expr {$dpsi2 + $b*$b*[lindex $::invM $i]}]
                }
    }

    # non-overlapping atoms (ids 4 and 16)
    foreach a [lindex $gradphi 0] b [lindex $gradpsi 3] {
                    set dphi2 [expr {$dphi2 + $a*$a*[lindex $::invM 0]}]
                    set dpsi2 [expr {$dpsi2 + $b*$b*[lindex $::invM 4]}]
                }

    return [list $dphi2 $dpsi2 $dphidpsi]
}

proc calc_difftensor_grad { args } {
    return {0.0 0.0 0.0}
}

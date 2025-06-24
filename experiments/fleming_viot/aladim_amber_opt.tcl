variable cos_coeff [list 11.833870900743202 -7.808474899351828 0.6449332357561219 -8.990607130973054 7.198940933336364 -6.95648845971968 -0.38621786071414765 -0.4433240832648866 -0.25483548766186925 -0.2290177491853248 0.14909564134174294 -0.047549855130133795 -0.14831926508217716 0.9467805938818169 -0.613769803799527]
variable sin_coeff [list 7.018759037141378 -26.021291153892335 4.337023790555449 -6.145919053914423 0.4052806810999208 2.144522922075821 -2.835129045048813 3.3731527691340855 -0.5587885599789834 0.49029182011932176 -0.7622692908687091 1.3197508962208049 -0.5866239215245678 0.6196352946467979 -0.637960388048086]


variable freqs [list 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 13.0 14.0 15.0]

proc calc_instate4 {args} {
    cv help
    set phi [cv colvar phi value]
    set psi [cv colvar psi value]

    set phim 53.0
    set psim 28.0

    # translate to closest periodic image

    set phi [expr {abs($phi - 360 - $phim) < abs($phi - $phim) ? $phi - 360 : $phi}]
    set psi [expr {abs($psi - 360 - $psim) < abs($psi - $psim) ? $psi - 360 : $psi}]

    set phi [expr {abs($phi + 360 - $phim) < abs($phi + $phim) ? $phi + 360 : $phi}]
    set psi [expr {abs($psi + 360 - $psim) < abs($psi - $psim) ? $psi + 360 : $psi}]

    # compute polar coordinates

    set r [expr {sqrt(($phi - $phim)*($phi - $phim) + ($psi - $psim)*($psi - $psim))}]
    set t [expr {atan2($psi - $psim, $phi - $phim)}]

    set R 132.37835147275467
    set i 0

    foreach f $::freqs {
        set R [expr {$R + [lindex $::cos_coeff $i]*cos($f*$t) + [lindex $::sin_coeff $i]*sin($f*$t)}]
        incr i
    }

    return [expr { $r < $R ? 1 : 0}]
}

proc calc_instate4_grad { args } { return 0.0 }
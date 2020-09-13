% ***********************************************************************
% calctav.m
% ***********************************************************************
% Stern F. (1964), Transmission of isotropic radiation across an
% interface between two dielectrics, Appl. Opt., 3(1):111-113.
% Allen W.A. (1973), Transmission of isotropic light across a
% dielectric surface in two and three dimensions, J. Opt. Soc. Am.,
% 63(6):664-666.
% ***********************************************************************
function tav=calctav(alfa,nr)

    rd  = pi/180;
    n2  = nr.^2;
    np  = n2+1;
    nm  = n2-1;
    a   = (nr+1).*(nr+1)/2;
    k   = -(n2-1).*(n2-1)/4;
    sa  = sin(alfa.*rd);

    b1  = (alfa~=90)*sqrt((sa.^2-np/2).*(sa.^2-np/2)+k);
    b2  = sa.^2-np/2;
    b   = b1-b2;
    b3  = b.^3;
    a3  = a.^3;
    ts  = (k.^2./(6*b3)+k./b-b/2)-(k.^2./(6*a3)+k./a-a/2);

    tp1 = -2*n2.*(b-a)./(np.^2);
    tp2 = -2*n2.*np.*log(b./a)./(nm.^2);
    tp3 = n2.*(1./b-1./a)/2;
    tp4 = 16*n2.^2.*(n2.^2+1).*log((2*np.*b-nm.^2)./(2*np.*a-nm.^2))./(np.^3.*nm.^2);
    tp5 = 16*n2.^3.*(1./(2*np.*b-nm.^2)-1./(2*np.*a-nm.^2))./(np.^3);
    tp  = tp1+tp2+tp3+tp4+tp5;
    tav = (ts+tp)./(2*sa.^2);

end

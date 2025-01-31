

% Kodun çalışabilmesi için InterX fonksiyonu kullanılmaktadır
% Stall etkisi gözetilerek
% grafik çizdirilmiştir.
clc;
clear;
alt=0:50:45000;
W=324720;
S=88.258;
CLmax=2.39;
for j = 1:length(alt)
    [temp, sos, press, rho, rhoh] = atmos_1976_SI(alt(j));
    V = linspace(30, 350, 1000);
    for i = 1:1000
        CL(i) = 2 * W / (rho * (V(i)^2) * S);
        M(i) = V(i) / sos;
        [CD0(i), K(i)] = drag_div_GulfstreamIV(M(i));
        CD(i) = CD0(i) + K(i) * CL(i)^2;
        TR(i) = 1/2 * rho * V(i)^2 * S * CD(i);
        [TA(i), ct(i)] = Two_Rolls_Royce_Tay_611(alt(j), V(i));
    end
    P = InterX([V;TR],[V;TA]); %InterX fonksiyonundan kesişim noktalarını çekmekte.
    if ~isempty(P)
        Vmin(j)=P(1,1);
        Vmax(j)=P(1,2);
        alt2(j)=alt(j);
        Vstall(j) = sqrt((2 * W) / (rho * S * CLmax));
        if Vstall(j)>Vmin(j)
            Vmin(j)=Vstall(j);
        else
           continue
        end
    else
       break
    end
end
plot(Vmin,alt2)
hold on
plot(Vmax,alt2)
plot([Vmin(end) Vmax(end)], [alt2(end) alt2(end)])
grid on
plot(Vmin, alt2, 'g', 'DisplayName', 'V_{min}');
plot(Vmax, alt2, 'b', 'DisplayName', 'V_{max}');
xlabel('Velocity (m/s)');
ylabel('Altitude (m)');
title('Flight Envelope');

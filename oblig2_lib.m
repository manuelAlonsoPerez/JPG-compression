function f = oblig2_lib
    
    f.DCT_2D=@DCT_2D;
    f.IDCT_2D=@IDCT_2D;
    f.my_entropy=@my_entropy;
    f.quantizier=@quantizier;
    f.dequantizier=@dequantizier;
    f.makehistograms= @makehistograms;

end

% Beregns vannlig, Normalisert og kumulativt histogramene
function [h_i,p_i,c_i] = makehistograms(img,greylevels)
    img = floor(img); 
    [M,N] = size(img);
    h_i = zeros(1,greylevels+1);
    c_i = zeros(1,greylevels+1);

    for i=1:M
       for j=1:N
           index = img(i,j)+1; % Maps 0 to index 1
           h_i(index)= h_i(index)+ 1;
       end
    end
    p_i=h_i/(N*M);

    c_i(1)= p_i(1);
    for i = 2 : greylevels
        c_i(i) = c_i(i-1) + p_i(i);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%          DCT blokkvis funksjonen               %%%
%%%% For hvert mxn blokk  beregnes DCT-en matrise,  %%%
%%%% multipliserer den punktvis med bilde blokken   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = DCT_2D(im,m,n)
   
    [M, N] = size(im);

    M8 = M + ( m-mod(M,m));
    N8 = N + ( n-mod(N,n));

    F = zeros(M8, N8);
    
    for u_offsett=1:m:M8-m
        for v_offsett=1:m:N8-n
            for u=1:m
                if (mod(u,m)==1) 
                    cu = 1/sqrt(2);       
                else
                    cu=1;  
                end
                for v=1:1+n
                    if (mod(v,n)==1) 
                        cv = 1/sqrt(2);  
                    else
                        cv=1;
                    end
                    D=zeros(m,n);
                    for x=0:7
                        for y=0:7
                             D(x+1,y+1) = cos(((u-1)*pi*(2*x+1))/16)*cos(((v-1)*pi*(2*y+1))/16);
                        end
                    end         
                    F(u+u_offsett-1,v+v_offsett-1)= sum(sum(D.*im(u_offsett:u_offsett+7,v_offsett:v_offsett+7)))*cu*cv/4;
                end
            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%          Kvantiserings funksjonen                %%%
%%%%  For hvert mxn blokk  dividerer DCT matrisen     %%%
%%%%  punktvis med  tilsvarende  bilde blokken        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function im = quantizier(im,Q)
    [M8, N8] = size(im);

    for u_offsett=1:8:M8-8
        for v_offsett=1:8:N8-8
            im(u_offsett:u_offsett+7,v_offsett:v_offsett+7)=(im(u_offsett:u_offsett+7,v_offsett:v_offsett+7)./Q);
        end
    end
    im=round(im);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%            Entropi  funksjonen                   %%%
%%%%  Beregnes den normaliserte histogrammet p(i)     %%%
%%%%  og entropien av verdiene der p(i)!=0            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ent = my_entropy(F)
    [M8, N8] = size(F);
    Max=max(max(F));
    Min=abs(min(min(F)));
    h=zeros(1,Max+Min+1);
    
    % Calculation of the histogram
    for i=1:M8
       for j=1:N8
           index = F(i,j);
           h(index+Min+1)= h(index+Min+1)+ 1;
       end
    end

    % Calculation of the Normalized histogram
    p=h/(M8*N8);

    % Calculation of the Entropy
    ent=0;
    for i = 1:Max+Min+1
        if(p(i)~=0)
            ent = ent - p(i)*log2(p(i));
        end
    end
    prob= sum(p);
    disp(prob)

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%       Invers Kvantiserings funksjonen                %%%
%%%%  For hvert mxn blokk  multipliserer DCT matrisen     %%%
%%%%   punktvis med  tilsvarende  bilde blokken           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function im = dequantizier(im,Q)
    [M8, N8] = size(im);

    for u_offsett=1:8:M8-8
        for v_offsett=1:8:N8-8
            im(u_offsett:u_offsett+7,v_offsett:v_offsett+7)=(im(u_offsett:u_offsett+7,v_offsett:v_offsett+7).*Q);
        end
    end
    im=round(im);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%          IDCT blokkvis funksjonen               %%%
%%%% For hvert mxn blokk  beregnes IDCT-en matrise,  %%%
%%%% multipliserer den punktvis med DCT blokken      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = IDCT_2D(F,m,n)
    [M8, N8] = size(F);
    f = zeros(M8, N8);

    for x_offsett=1:m:M8-8
        for y_offsett=1:n:N8-8

            for x=1:m
                for y=1:1+n

                    D=zeros(m,n);
                    for u=0:7
                        for v=0:7
                            if (u==0) 
                                cu = 1/sqrt(2);       
                            else
                                cu=1;  
                            end
                            if (v==0) 
                                cv = 1/sqrt(2);  
                            else
                                cv=1;
                            end
                            D(u+1,v+1) = cu*cv* cos( (u*pi*(2*(x-1)+1))/16 )...
                                * cos( (v*pi*(2*(y-1)+1))/16 );
                        end
                    end         
                    f(x+x_offsett-1,y+y_offsett-1)= sum(sum(F(x_offsett:x_offsett+7,y_offsett:y_offsett+7).*D))/4;
                end
            end
        end
    end
    

end


package falconsTools;
import java.util.*;
public class TailDesign {	
	static double Sw, Sf, Sh, Sref ; 
	static double aw, af, ah, awf ,Eslope, Cmslope ; 
	static double x_cg, x_ac_f, x_ac_w, x_np, x_ac_h, fuselageArm, wingArm , HstabArm ; 
	static double Cw, Cf,  Ch , Cref , bh ; 
	static double CL0w, CL0f ;
	static double AOIf, AOIw, AOIh, E, E0;
	static double HstabEfficiency = 0.9  ;
	static double Cmf, Cmw, Cmh, Cmref ,Cm0ref ;
	static double AoAw , AoAf , AoAh , AoAwf, AoAref ;
	static double ARwf , pi = 3.14159265358 , density = 1.225 , mu = 1.86*Math.pow(10 ,  5 ) ;
	static double CLw , CLf , CLh, CLwf , CLref ;
	static double fuselageVolume , wingVolume, HstabVolume ; 	
	static double ah_2d = 0.15 , maxthickness = 12 , maxthicknesspos = 29 ; 
	
	static double Cv, bv, lambdav, offsetv, VstabArm ;
	static double Sv , MACv , ARv, ev , av , CLv , VstabVolume ;
	static double bw, bf, bwf ; 
	static double CnB , CnBwf , Sigmaslope;
	
	//AR 2.9 - 3.3
	//optimum tail moment arm 1 - 1.5
	//static margin 0.1 - 0.15
	
	public static void input () { 
		
		//INPUT
		Sw = 1.41 ; 
		Sf = 0.33 ; 
		Cw = 0.56 ; 
		bw = 0 ; 
		bf = 0 ;
		lambdav = 0 ; 
		offsetv = 0 ;
		VstabArm = 0 ;
		Cf = 0.831 ; 
		aw = 0.07114 ; 
		af = 0.05776  ; 
		
		Sref = Sw + Sf; 
		
		x_cg = 0.375; 
		x_ac_f = 0.235 ;
		x_ac_w = 0.375 ;
		
		Cref = ( Cw*Sw + Cf*Sf  ) / Sref ; 
		awf = ( aw*Sw + af*Sf ) / Sref ;
		ARwf = 3 * 3 / Sref ; 
		CL0w = 1.01 ;
		CL0f = 0.39 ;
		AOIf = 3 ; 
		AOIw = 5 ;
		
		fuselageArm = x_cg - x_ac_f ;
		wingArm = x_cg - x_ac_w ; 
		Cmf = -0.04 ;
		Cmw = -0.068 ;
		//Cmf = -0.1091181818 ; 
		//Cmw = -0.13858 ;		
	}	
	
	public static void main ( String args[] ) { 
		input() ; 
				
		fuselageVolume = ( Sf ) / ( Sref * Cref ) ;
		wingVolume = ( Sw  ) / ( Sref * Cref ) ;
		
		for ( HstabArm = -1 ; HstabArm >= -1.3 ; HstabArm -= 0.05 ) { 
			for ( Ch = 0.3 ; Ch <= 0.45 ; Ch += 0.005 ) { 
				for ( bh = 0.5 ; bh <= 1.3 ; bh += 0.005 ) { 
					for ( AOIh = -2 ; AOIh >= -5 ; AOIh -= 0.5 ) { 
						
						double ARh = bh / Ch ;
						if (ARh < 2.8 ) 
							continue ; 						
						double eh = calculate_e_Cd0("e", Ch, bh, AOIh, 1, 1, 13) ; 
						ah = ah_2d / ( 1 + 57.3 * ah_2d / ( pi * eh * ARh)) ;				
						AoAw = AOIw ; 
						AoAf =  AOIf ; 
						AoAwf = ( AoAw*Sw + AoAf*Sf ) / Sref ;
						//E0 = 2* CLwf / ( pi * ARwf) ; 
						//Eslope = 2 * awf / ( pi * ARwf) ; 
						E = E0 + Eslope*AoAwf ;
						AoAh =  AOIh - E ;
						Sh = bh * Ch ;

						HstabVolume =  ( Sh * HstabArm ) / ( Sref * Cref ) ;
						CLf = CL0f + af*AoAf ;
						CLw = CL0w + aw*AoAw ;
						CLh = ( ah * (1 - Eslope)*(AoAh)) ;
						CLwf =( CLw*Sw + CLf*Sf ) / Sref ; 
					
						//CM calculation
						
						double fuselageContribution = fuselageVolume*Cf*Cmf + fuselageVolume*fuselageArm*CLf ;
						double wingContribution = wingVolume*Cw*Cmw + wingVolume*wingArm*CLw ;
						double tailContribution = HstabEfficiency*HstabVolume*CLh ;								
												
						Cm0ref = fuselageContribution + wingContribution + tailContribution ; 
						
						//CM slope calculation 
												
						fuselageContribution = fuselageVolume*fuselageArm*af ;
						wingContribution = wingVolume*wingArm*aw;
						tailContribution = HstabEfficiency*HstabVolume*( 1 - Eslope)*ah ;
						Cmslope = fuselageContribution + wingContribution + tailContribution ; 	
												
						double trimAoA = -Cm0ref / ( Cmslope ) ;		
						double staticMargin =- (fuselageContribution + tailContribution) / ( ( aw * Sref )/Sw ) ;	
						double neutralPoint = x_ac_w + staticMargin* Cref ; 		
						
						//Vstab
						/*for ( bv = 0.2; bv <= 0.5 ; bv += 0.02 ) { 
							for ( lambdav = 0.2 ; lambdav <= 1 ; lambdav += 0.02 ) {
								Cv = Ch ;
								Sv = Cv * ( 1 + lambdav ) / 2 * bv ;
								MACv = Cv * 2/3 * (( 1 + lambdav + lambdav*lambdav ) ÷ ( 1 + lambdav )) ;
								VstabArm = HstabArm - 0.25 * Ch + Cv - 0.75 * MACv ;
								VstabVolume = VstabArm * Sv / ( Swf * bwf ) ;
								
								CnBv = CnBwf + efficiencyv * VstabVolume * av * ( 1 + Sigmaslope ) ;
								
							}
						} */
					
						
						if ( staticMargin > 0.1 && trimAoA > 0 && trimAoA < 1  && HstabVolume > 0.2 ) {
							System.out.println("Vh " + HstabVolume ) ;		
							System.out.println("c " + Ch ) ;		
							System.out.println("b " + bh ) ;		
							System.out.println(  "ARh" + ARh );
							System.out.println("AOI " + AOIh ) ;		
							System.out.println("tail arm " + HstabArm) ;		
							System.out.println("static Margin " + staticMargin*100) ;	
							System.out.println("Cm0 " + Cm0ref) ;
							System.out.println("Cmslope" + Cmslope) ;
							System.out.println ( "trim AoA :" + trimAoA) ;
							System.out.println ( "E0:" + E0 ) ;
							System.out.println(  (0.235 + 0.14 - HstabArm - 0.25 * Ch ) );
							System.out.println() ; 
						}
					}
				}
			}
		}		
	}
	public static double calculate_e_Cd0(String choice , double c , double b , double aoi , double phi , double lambda , double Velocity ) {
		double Ne = 0 ; 
		double S1 = b * c * phi ;
		double Ma = 0.04 ;
		double S2 = b * c * ( 1 + lambda ) /2 * ( 1 - phi ) ;
		double Sref = S1 + S2 ; 
		double ceff = Sref / b ;
		double Re = ceff * density * Velocity / mu ; 
		double AR = b * b / Sref ; 
		double TR = ( S1 + S2 * lambda ) / Sref ; 
		double flambda =  0.005* ( 1 + 1.5 * Math.pow ( (TR - 0.6) , 2 ) ) ;
		double sweepAngle = Math.atan((-0.25 * c + 0.25 * lambda * c) / (b / 2));
		double e = 1/((1 + 0.12 * Ma* Ma) *
                ((1 + (0.142 + flambda * AR * Math.pow(10 * maxthickness / 100, 0.33)) / Math.pow(Math.cos(sweepAngle), 2) +
                (0.3 * Ne + 0.1) / Math.pow(4 + AR, 0.8)))) ;		
		double Cf = 0.074 / Math.pow( Re,0.2 ) ;
		double FF = ((1 + ( maxthicknesspos / 100 > 0.3 ? 1.2 : 2) * maxthickness / 100 + 100 * Math.pow(maxthickness / 100, 4)) * 1.05);
		double Swet = 2*(1+0.25*maxthickness/100) ;
		double Cd0 = Swet * FF * Cf ;
		if ( choice.equals ( "e" ) ) 
			return e ; 
		else return Cd0 ;
	}
}

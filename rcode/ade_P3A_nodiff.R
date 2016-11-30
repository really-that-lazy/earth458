rm(list=ls())
graphics.off()
##############################################################
#Written by Jason Davison
#C is the contaminant or paramater being solved in the equaion
##############################################################

##############################################################
#Model Parameters
##############################################################
Cr = 1.0	   #Courant Number
L = 10		   #Length of Problem (m) Make it a little longer than the point of interest
dz = ceiling(L/10)/10  #Discretization (m)
R = 3.6             #Retardation Factor
##disp = 8.2*10^-10/R      #Dispersion constant (m^2/s)
disp = 0
u = 0 /R	   #Velocity Left to Right (m/s) 
Runtime = 5*10^10 #1.5*L/u   #Solution Durration in (s) May need to edit this value to get correct values
##############################################################
#Boundary Conditions
##############################################################
###---Left Hand Conditions---#################################
LeftHandBC = 1  #Left hand boundary conditions
		#1 = Constant Concentration
		#2 = Free Flux Boundary (i.e. open boundary)
LeftHandC = 100   #Left hand boundary concentration (only for constant concentration conditions)
###---Right Hand Conditions---################################
RightHandBC = 2  #Right hand boundary conditions
                 #1 = Constant Concentration
                 #2 = Free Flux Boundary (i.e. open boundary)
RightHandC = 0   #Right hand boundary condition (only for constant concentration conditions)

##############################################################
#Time Stepping Parameters
##############################################################
#dt = (Cr*dz/u)/20  #Time step
dt = Runtime/800
inter = Runtime/dt      #Number of simultation timesteps

##############################################################
#Initial Conditions
##############################################################
BackgroundC = 0  #Background Concentration for simulation (mg/L)

StartOfPlume = 1 #This is the start of the Plume (m) 
EndOfPlume = 2    #This is the ending value of the Plume (m) 
PlumeC = 0       #Concentration of Plume (mg/L)
ConstantPlume = 2   #Constant Source Concentration
		    #1 = Plume is Constant Source (i.e. it will always be there)
		    #2 = Plume is Instaneous Source

C = matrix(0,(L/dz),inter+1)	
if(PlumeC>0){
  C[1:(L/dz),1] = rep(BackgroundC,L/dz)
  C[floor(StartOfPlume/dz):floor(EndOfPlume/dz),1] = rep(PlumeC,1+floor(EndOfPlume/dz)-floor(StartOfPlume/dz))
}
##############################################################
#Plotting Output
##############################################################
simpleplot = TRUE  #Plots all time values
prettyplot = TRUE   #Plots specific values

Plot_Time = c(Runtime/10,3*Runtime/10,Runtime/2,7*Runtime/10, 9*Runtime/10)  #Output times for pretty plots must be less than Runtime
#Plot_Time = c(10^3, 2*10^3)
#Plot_Location = c(L/5,2*L/5,3*L/5,4*L/5)          #Must be greater than 0
Plot_Location = c(1,10)
##################
if(simpleplot){X11()}
##################
#Time Loop
##################
for(i in 1:inter){
	######-----Do Not Edit Solver----#####
	a = -dt/dz^2*disp-u*dt/(dz)
	b = 1+2*dt/dz^2*disp+u*dt/(dz)
	c = -dt/dz^2*disp
	###############
	# A Matrix
	###############
	A_a = diag(a,L/dz-2)
	A_b = diag(b,L/dz)
	A_c = diag(c,L/dz-2)

	A = A_b
	A[2:(L/dz-1),1:(L/dz-2)]  = A[2:(L/dz-1),1:(L/dz-2)] + A_a
	A[2:(L/dz-1),3:(L/dz)]  = A[2:(L/dz-1),3:(L/dz)] + A_c
	A[1,1] = 1
	A[L/dz,L/dz] =1
	
	###################
	#Boundary Conditions
	###################
	if(LeftHandBC==1){C[1,i] = LeftHandC}
	if(LeftHandBC==2){C[1,i] = C[2,i]}
        if(RightHandBC==1){C[L/dz,i] = RightHandC}
	if(RightHandBC==2){C[L/dz,i] = C[L/dz-1,i]}	
        if(ConstantPlume==1){C[(StartOfPlume/dz):(EndOfPlume/dz),i] = rep(PlumeC,1+EndOfPlume/dz-StartOfPlume/dz)}

	#B vector for back solve solution
	B = C[,i]
	C[,1+i] = solve(A,B)	#New solution

        ######-----Do Not Edit Solver----#####
	#simple plot:
        x_length = c(1:(length(C[,i])))*dz
	if(simpleplot){plot(x_length,C[,i+1], type="o", col="blue",xlim=c(0,L),ylim=c(0,max(C[,1])+0.2), xlab="Length (m)", ylab="Concentration (mg/L)")}
}

##########################
#----Ploting Utility----##
##########################
if(prettyplot){
number_output = length(Plot_Time)
graphics.off()
xrange <- range(x_length)
yrange <- range(0,(max(C[,1])*1.25))

pdf(file="ADE_CDprofile_a4q3.pdf")
par(mar=c(6,6,2,2))
plot(xrange, yrange, type="n", ylab = " ", xlab = " ", main = " ", cex.axis = 2)
mtext("Distance (m)", side = 1, line = 2.5, cex = 2)
mtext(~"Concentration (mg/L)", side = 2, line = 2.5, cex = 2)
colors <- topo.colors(number_output)
linetype <- c(1:3)
plotchar <- c(1:number_output)

for (i in 1:number_output) { 
plotchar[i] = toString(format(signif(Plot_Time[i],digits=2),scientific=TRUE))
  lines(x_length, C[,floor((Plot_Time[i]/dt))], type="l", lwd=3, lty=linetype[1], col=colors[i], pch=plotchar[1])
}
# add a legend
legend('topright', plotchar, cex=1.0, col=colors, lty=linetype[1], lwd = 3, bg = "white",ncol=3)
dev.off()

#For breakthrough curves
time_vector = c(0:(Runtime/dt))*dt
number_output = length(Plot_Location)
graphics.off()
xrange <- range(0,Runtime)
yrange <- range(0,(max(C[,1])*1.25))

pdf(file="ADE_BTcurve_a4q3_no.pdf")
par(mar=c(6,6,2,2))
plot(xrange, yrange, type="n", ylab = " ", xlab = " ", main = " ", cex.axis = 2)
mtext("Time (s)", side = 1, line = 2.5, cex = 2)
mtext(~"Concentration (mg/L)", side = 2, line = 2.5, cex = 2)
colors <- topo.colors(number_output)
linetype <- c(1:3)
plotchar <- c(1:number_output)

for (i in 1:number_output) {
  plotchar[i] = toString(signif(Plot_Location[i],digits=3))
  lines(time_vector, C[floor((Plot_Location[i]/dz)),], type="l", lwd=3, lty=linetype[1], col=colors[i], pch=plotchar[1])
}
# add a legend
legend('topright', plotchar, cex=1., col=colors, lty=linetype[1], lwd = 3, bg = "white",ncol=3)
dev.off()
}

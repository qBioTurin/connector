#' @title Run docker container
#' @description This is an internal function executing a docker container.
#' @return 0 if success, 1 if dockerid file is present, 2 if docker execution fails.
#' @author Beccuti Marco, Pernice Simone
#'
#' @import utils
#' @examples
#'\dontrun{
##'     #running runDocker
#'      docker.run(params=NULL)
#'
#' }
#' @export
docker.run <- function(){

	# to check the Docker ID by file
	if (file.exists("dockerID")){
		cat("\n\nDocker does not start, there is already a docker container running che dockerID file!!!\n\n")
		system("echo 1 > ExitStatusFile 2>&1")
		return(1)
	}


	## to execute docker
	cmdl = paste("docker run --cidfile=dockerID --privileged=true -e DISABLE_AUTH=true -p 8787:8787 -d qbioturin/connector:latest", "\n\n", sep="")
	cat(cmdl)
	system(cmdl)

	## Get the Docker ID from file
	dockerid=readLines("dockerID", warn = FALSE)
	cat("\nDocker ID is:\n",substr(dockerid,1, 12),"\n")

	browseURL(url = "http://localhost:8787/")

	## Check the Docker container status
	dockerStatus=system(paste("docker inspect -f {{.State.Running}}",dockerid),intern= T)
	while(dockerStatus=="true"){
		Sys.sleep(10);
		dockerStatus=system(paste("docker inspect -f {{.State.Running}}",dockerid),intern= T)
		cat(".")
	}
	cat(".\n\n")
	## Check the Docker container exit status
	dockerExit <- system(paste("docker inspect -f {{.State.ExitCode}}",dockerid),intern= T)
	cat("\nDocker exit status:",dockerExit,"\n\n")

	if(as.numeric(dockerExit)!=0){
		system(paste("docker logs ", substr(dockerid,1,12), " &> ", substr(dockerid,1,12),"_error.log", sep=""))
		cat(paste("\nDocker container ", substr(dockerid,1,12), " had exit different from 0\n", sep=""))
		cat("\nExecution is interrupted\n")
		cat("The container's log is saved at: ")
		system(paste0("docker inspect --format=","'{{.LogPath}}' ",dockerid))
		cat(paste("Please send to beccuti@unito.it this error: Docker failed exit 0,\n the description of the function you were using and the following error log file,\n which is saved in your working folder:\n", substr(dockerid,1,12),"_error.log\n", sep=""))
		system("echo 2 > ExitStatusFile 2>&1")
		return(3)
	}
	else
	{
		file.remove("dockerID")
		system(paste("docker rm -f ",dockerid),intern= T)
	}

	#Normal Docker execution
	system("echo 0 > ExitStatusFile 2>&1")
	return(0)
}

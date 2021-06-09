package org.reactome.immport.ws.service;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.ServerSocket;
import java.net.Socket;

import javax.annotation.PreDestroy;

import org.apache.log4j.Logger;

public class RScriptService {
    private static Logger logger = Logger.getLogger(RScriptService.class);

    private String dataDir;
    private String script;
    private int port;
    private String workDir;
    private String command;
    // To kill it
    private Process process;

    public RScriptService() {
    }

    public int getPort() {
        return this.port;
    }

    public String getCommand() {
        return command;
    }

    public void setCommand(String command) {
        this.command = command;
    }

    public String getDataDir() {
        return dataDir;
    }

    public void setDataDir(String dataDir) {
        this.dataDir = dataDir;
    }

    public String getScript() {
        return script;
    }

    public void setScript(String script) {
        this.script = script;
    }

    public String getWorkDir() {
        return workDir;
    }

    public void setWorkDir(String workDir) {
        this.workDir = workDir;
    }
    
    @PreDestroy
    public void stopService() {
        if (process != null && process.isAlive())
            process.destroy();
    }

    /**
     * Start the R REST service by using Java process calling
     */
    public void startService() {
        try {
            // Get an available port
            ServerSocket socket = new ServerSocket(0);
            this.port = socket.getLocalPort();
            socket.close();
            // Build an asynchronous call
            String[] parameters = {
                    command,
                    workDir + File.separator + script,
                    port + "",
                    dataDir,
                    workDir
            };
            logger.info("Starting the R service...");
            logger.info("Working dir: " + workDir);
            ProcessBuilder builder = new ProcessBuilder(parameters);
            builder.directory(new File(workDir));
            builder.redirectErrorStream(true);
            // For the time being, this will go to system.out.
            // TODO: redirect to logger as info.
            //            builder.redirectOutput(Redirect.INHERIT);
            process = builder.start();
            // Merge the output into logger
            Thread t = new Thread() {
                public void run() {
                    try {
                        BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
                        String line = null;
                        while ((line = reader.readLine()) != null) {
                            logger.info(line);
                        }
                    }
                    catch(IOException e) {
                        logger.error(e.getMessage(), e);
                    }
                }
            };
            t.setPriority(t.getPriority() - 1);
            t.start();
            // Make sure the service start 
            // Poll the port to make sure it is reachable
            long time1 = System.currentTimeMillis();
            while (true) {
                try {
                    Thread.sleep(1000);
                    // Try again
                    Socket testSocket = new Socket("localhost", port);
                    if (testSocket.isConnected()) {
                        testSocket.close();
                        logger.info("The R service has started.");
                        break;
                    }
                    else
                        testSocket.close();
                }
                catch(Exception e) { // Do nothing
                    // Avoid to throw ConnectException, which should be thrown during starting up.
//                    logger.error(e.getMessage(), e);
                }
                long time2 = System.currentTimeMillis();
                // It is observed that the service starts slow when it starts for the first time after download under Mac.
                if ((time2 - time1) > 10000) { // Give 10 seconds for starting the R service
                    logger.info("The R service cannot start in 10 seconds. Please check!");
                    break;
                }
            }
        }
        catch(Exception e) {
            logger.error(e.getMessage(), e);
        }
    }
}

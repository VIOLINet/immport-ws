package org.reactome.immport.ws.service;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import org.apache.log4j.Logger;

public class RScriptService {
    private static Logger logger = Logger.getLogger(RScriptService.class);
    
    private String dataDir;
    private String script;
    private String port;
    private String workDir;
    private String command;
    
    public RScriptService() {
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

    public String getPort() {
        return port;
    }

    public void setPort(String port) {
        this.port = port;
    }

    public String getWorkDir() {
        return workDir;
    }

    public void setWorkDir(String workDir) {
        this.workDir = workDir;
    }
    
    /**
     * Start the R REST service by using Java process calling
     */
    public void startService() {
        String[] parameters = {
                command,
                workDir + File.separator + script,
                port,
                dataDir,
                workDir
        };
        Thread t = new Thread() {
            public void run() {
                try {
                    logger.info("Starting R service...");
                    Runtime runtime = Runtime.getRuntime();
                    Process process = runtime.exec(parameters);
                    InputStream is = process.getInputStream();
                    String output = output(is);
                    logger.info(output);
                    logger.info("R service started: " + process);
                }
                catch(Exception e) {
                    logger.error(e.getMessage(), e);
                }
            }
        };
        t.start();
    }
    
    private String output(InputStream is) throws IOException {
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line = null;
        StringBuilder builder = new StringBuilder();
        while ((line = br.readLine()) != null) {
            builder.append(line).append("\n");
        }
        br.close();
        isr.close();
        return builder.toString();
    }

}

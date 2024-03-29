package org.reactome.immport.ws.config;

import static org.hibernate.cfg.AvailableSettings.C3P0_ACQUIRE_INCREMENT;
import static org.hibernate.cfg.AvailableSettings.C3P0_MAX_SIZE;
import static org.hibernate.cfg.AvailableSettings.C3P0_MAX_STATEMENTS;
import static org.hibernate.cfg.AvailableSettings.C3P0_MIN_SIZE;
import static org.hibernate.cfg.AvailableSettings.C3P0_TIMEOUT;
import static org.hibernate.cfg.AvailableSettings.DRIVER;
import static org.hibernate.cfg.AvailableSettings.HBM2DDL_AUTO;
import static org.hibernate.cfg.AvailableSettings.PASS;
import static org.hibernate.cfg.AvailableSettings.SHOW_SQL;
import static org.hibernate.cfg.AvailableSettings.URL;
import static org.hibernate.cfg.AvailableSettings.USER;

import java.util.Properties;

import org.reactome.immport.ws.service.ImmportServiceConfig;
import org.reactome.immport.ws.service.RScriptService;
import org.reactome.immport.ws.service.ReactomeAnalysisConfig;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.ComponentScan;
import org.springframework.context.annotation.ComponentScans;
import org.springframework.context.annotation.Configuration;
import org.springframework.context.annotation.PropertySource;
import org.springframework.core.env.Environment;
import org.springframework.orm.hibernate5.HibernateTransactionManager;
import org.springframework.orm.hibernate5.LocalSessionFactoryBean;
import org.springframework.transaction.annotation.EnableTransactionManagement;

@Configuration
@PropertySource("classpath:application.properties")
@EnableTransactionManagement // Need this level of transaction manager. Otherwise, DAO cannot get transaction.
@ComponentScans(value = { @ComponentScan("org.reactome.immport.ws.service")})
public class AppConfig {

	@Autowired
	private Environment env;

	@Bean
	public ReactomeAnalysisConfig getReactomeAnalysisConfig() {
		ReactomeAnalysisConfig config = new ReactomeAnalysisConfig();
		config.setReactomeAnalysisURL(env.getProperty("reactome.analysis.service"));
		config.setReactomeFIServiceURL(env.getProperty("reactome.fi.service"));
		config.setReactomePathwayHierarchyURL(env.getProperty("reactome.pathway.hierarchy.service"));
		config.setModuleColors(env.getProperty("reactome.fi.moduleColors").split(";"));
		return config;
	}

	@Bean
	public ImmportServiceConfig getImmportServiceConfig() {
		ImmportServiceConfig config = new ImmportServiceConfig();
		config.setBiosampleMetatdataFileLocation(env.getProperty("biosample.metatdata.file"));
		config.setTestDiffGeneExpFileLocation(env.getProperty("test.diff.exp.result.file"));
		return config;
	}

	@Bean
	public LocalSessionFactoryBean getSessionFactory() {
		if (env.getProperty("is.deployed").equals("true")) return null;
		LocalSessionFactoryBean factoryBean = new LocalSessionFactoryBean();

		Properties props = new Properties();
		// Setting JDBC properties
		props.put(DRIVER, env.getProperty("mysql.driver"));
		props.put(URL, env.getProperty("mysql.url"));
		props.put(USER, env.getProperty("mysql.user"));
		props.put(PASS, env.getProperty("mysql.password"));
		// The following code prevents JDBC connection.
		//      props.put("mysql.requireSSL", env.getProperty("my.requireSSL"));
		//      props.put("mysql.verifyServerCertificate", env.getProperty("mysql.verifyServerCertificate"));
		//      props.put("mysql.useSSL", env.getProperty("mysql.useSSL"));

		// Setting Hibernate properties
		props.put(SHOW_SQL, env.getProperty("hibernate.show_sql"));
		props.put(HBM2DDL_AUTO, env.getProperty("hibernate.hbm2ddl.auto"));
		// The following cannot work. Have to set in the running property as the following:
		// -Dhibernate.dialect.storage_engine=innodb
		//      props.put(STORAGE_ENGINE, env.getProperty("hibernate.dialect.storage_engine"));

		// Setting C3P0 properties
		props.put(C3P0_MIN_SIZE, env.getProperty("hibernate.c3p0.min_size"));
		props.put(C3P0_MAX_SIZE, env.getProperty("hibernate.c3p0.max_size"));
		props.put(C3P0_ACQUIRE_INCREMENT, 
				env.getProperty("hibernate.c3p0.acquire_increment"));
		props.put(C3P0_TIMEOUT, env.getProperty("hibernate.c3p0.timeout"));
		props.put(C3P0_MAX_STATEMENTS, env.getProperty("hibernate.c3p0.max_statements"));

		factoryBean.setHibernateProperties(props);
		factoryBean.setPackagesToScan("org.reactome.immport.ws.model");

		return factoryBean;
	}

	@Bean
	public HibernateTransactionManager getTransactionManager() {
		if (env.getProperty("is.deployed").equals("true")) return null;
		HibernateTransactionManager transactionManager = new HibernateTransactionManager();
		transactionManager.setSessionFactory(getSessionFactory().getObject());
		return transactionManager;
	}

	@Bean
	public RScriptService getRService() {
		RScriptService service = new RScriptService();
		service.setDataDir(env.getProperty("r.data.dir"));
		service.setScript(env.getProperty("r.script"));
		service.setWorkDir(env.getProperty("r.work.dir"));
		service.setCommand(env.getProperty("r.command"));
		service.startService();
		return service;
	}
}
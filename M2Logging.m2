newPackage(
     "M2Logging",
     Version => "0.1", 
     Date => "30.07.2013",
     Authors => {{
           Name => "Jakob Kroeker", 
           Email => "kroeker@uni-math.gwdg.de", 
           HomePage => "http://www.crcg.de/wiki/User:Kroeker"}  
      },
     Configuration => {},
     Headline => "a rudimentary logger implementation (not threadsafe)",
     DebuggingMode => true
)

needsPackage "SimpleDoc";
needsPackage "Text";

-- todo: have Problems to export a symbol ... 




-- Merke: um variablen vor Auswertung zu schützen, setze diese in ""; dies ist z.B. für die Dokumentation von "LogLevel" notwendig!

export {
    LogLevel, 
    Logger,
    "M2LoggerExport"
   --"getLogLevel"
}

LogLevel = new HashTable from  {
FATAL   => 5,
ERROR   => 4,
WARNING => 3,
INFO    => 2,
DEBUG   => 1,
TRACE   => 0,
};

M2LoggerProtect = ()->
(
    protect setLogLevel;
    protect getLogLevel;
    protect getLogLevelName;
    protect name;
    protect displayName;
    protect displayMsgLevel;

    protect fatal;
    protect warning;

    protect FATAL;
    protect ERROR;
    protect WARNING;
    protect INFO;
    protect DEBUG;
    protect TRACE;
);


M2LoggerExport = ()->
(
   exportMutable( setLogLevel );
    exportMutable( getLogLevel );


   --getLogLevel1 := global getLogLevel;
   --getLogLevel1 = getLogLevel;
   --getLogLevel1 = getGlobalSymbol(M2Logging.Dictionary, "getLogLevel");

   exportMutable( getLogLevelName );
   exportMutable( name );
   exportMutable( displayName );
   exportMutable( displayMsgLevel );

   exportMutable( fatal );
   exportMutable( warning );
 
   exportMutable( FATAL );
   exportMutable( ERROR );
   exportMutable( WARNING );
   exportMutable( INFO );
   exportMutable( DEBUG );
   exportMutable( TRACE );
);

getLogLevel := getLogLevel;

M2LoggerExport(); 





--eventually create a hash table with existing loggers (hash by name?)

loggerTable := new MutableHashTable;

 
LogLevelName := (new  HashTable from apply(keys LogLevel, key -> LogLevel#key=>toString key));




-- returns (or creates) a logger which is identified by name in a HashTable
  
LogLevelIsValid := method();

LogLevelIsValid( ZZ ) := Boolean => (  level )->
(
   for key  in keys LogLevel do
   (
     if level==LogLevel#key then return true;
   );
   return false;
);


-- todo : is not threadsafe.

Logger = method();

Logger(String, ZZ) := HashTable => ( loggerName, pLevel )->
(
   if not LogLevelIsValid(pLevel) then 
         error("invalid log level; see valid levels in 'LogLevel'") ;

    logger :=  null;
    if loggerTable#?loggerName then 
    (
        logger = loggerTable#loggerName;
        if  not logger.getLogLevel()== pLevel then 
            print ("warning: updating log level of logger " |  loggerName );
        logger.setLogLevel( pLevel );
        return logger;
    );
    --
    logger = new MutableHashTable;
    --
  
  
    level := null;
    logger.name = loggerName;
    logger.setLogLevel = method();
    --
    logger.setLogLevel(ZZ) := ( plevel )->
    ( 
       if not LogLevelIsValid(plevel) then 
           error("invalid log level; see valid levels in 'LogLevel'") ;
       level = plevel;
    );

    logger.displayName = method();
    --
    mDisplayName := true;
    logger.displayName(Boolean) := ( bDisplayName )->
    ( 
      mDisplayName = bDisplayName;
    );

    mDisplayMsgLevel := true ;
    logger.displayMsgLevel = method();

    logger.displayMsgLevel(Boolean) := ( bDisplayMsgLevel )->
    ( 
      mDisplayMsgLevel = bDisplayMsgLevel;
    );


    logger.getLogLevel  = ()-> return level;
    logger.getLogLevelName  = ()-> return LogLevelName#level;
 
    --    
    logger.log = method();
    logger.log (ZZ,String) := (msgLvl,msg )->
    (
       if (msgLvl>=level) then
       (
           dspmsg := "";
           if mDisplayName then dspmsg= dspmsg |  logger.name  | ":"; 
           if mDisplayMsgLevel then dspmsg= dspmsg | LogLevelName#level  | ":";

           dspmsg  =  dspmsg | " " | msg ;
           print( "-- " | dspmsg   );
       );
    );

    logger.fatal = method();
    logger.fatal(String) := ( msg )-> ( logger.log(LogLevel.FATAL, msg ); );
    logger.error = method();
    logger.error (String) := ( msg )-> ( logger.log(LogLevel.ERROR, msg ); );
    logger.warning = method();
    logger.warning (String) := ( msg )-> ( logger.log(LogLevel.WARNING, msg ); );
    logger.info = method();
    logger.info (String) := ( msg )-> ( logger.log(LogLevel.INFO, msg ); );
    logger.debug = method();
    logger.debug (String) := ( msg )-> ( logger.log(LogLevel.DEBUG, msg ); );
    logger.trace = method();
    logger.trace (String) := ( msg )-> ( logger.log(LogLevel.TRACE, msg ); );


    --
    constLogger := new HashTable from logger;

    logger.setLogLevel( pLevel );
    
    loggerTable#loggerName= constLogger;
    return constLogger;
);

Logger(String) := HashTable => ( loggerName )->
(
    if loggerTable#?loggerName then 
    (
        return loggerTable#loggerName;
    );
    return Logger(loggerName, LogLevel.ERROR );
);

--
doc ///
    Key
       Logger
    Headline
       get a non-threadsafe logger 
    Description
        Text
             A Logger is constructed (or obtained by name )
        Example
            ll = Logger("exampleLogger");
            ll.setLogLevel(LogLevel.WARNING);
            ll.warning("this is a warning");
            ll.info("the info is not displayed, since the logLevel is too restrictive");
            ll.setLogLevel( LogLevel.INFO );
            ll.info("now the info is  displayed ");
        Text
            ll.Log( LogLevel.DEBUG )
///


doc ///
    Key
       "LogLevel"
    Headline
          List of valid log levels
    Description
        Text
            represents a list of valid log levels 
        Example
            LogLevel
///




doc ///
    Key
        M2Logging
    Headline
          A non-threadsafe logger package
    Description
        Text
            Implements a non-threadsafe logger package \break
            \break
            A Logger is constructed (or obtained by name via )
            \,\,  \bullet \,   @TO Logger@ ({\tt (loggerName, "loglevel" ) })  \break     
///


-- restart
-- loadPackage "M2Logging"
TEST /// 
    debug M2Logging

    M2LoggerProtect()

    logger = Logger( "testLogger",5 );
    assert( logger.getLogLevel()==5 );
    logger.setLogLevel( LogLevel.INFO );
    assert( logger.getLogLevel()==LogLevel.INFO );
    logger.warning("test");
    logger.displayName(true);
    logger.displayMsgLevel(true);
    logger.log( LogLevel.INFO, " testmsg ");

///

end


function UTLsbnserver (portnum)
  
  QUITMESSAGE    = "quit UTLsbnserver";
  CONFIRMMESSAGE = "confirmed";

  ## CREATE THE SOCKET AND WAIT FOR CONNECTIONS
  s = socket(AF_INET, SOCK_STREAM, 0);
  if s < 0
    error("cannot create socket\n");
  end
  
  if bind(s, portnum) < 0
    error("bind failed\n");
  end

  if listen(s, 1) < 0
    error("listen failed\n");
  end

  ##MAIN LOOP
  while 1

    ##ACCEPT CONNETCIONS
    c = accept(s);
    if c < 0
      error("connection error")
    end
    
    ## READ COMMANDS FROM THE SOCKET  
    msg = readstring (c)
    
    ##IF CLIENT SENT SHUTDOWN MESSAGE EXIT
    if strcmp (msg,QUITMESSAGE)
      printf("client requested server shutdown, goodbye!\n");
      disconnect(c); disconnect(s);
      break
    end
    
    ##EXECUTE COMMANDS FROM THE CLIENT
    [A,B,C] = eval(msg);
    
    ##SEND OUTPUT TO THE CLIENT
    str = [ sprintf("%17g ",A) "\n" sprintf("%17g ",B)...
	   "\n" sprintf("%17g ",C) "\n"]
    
    send(c,str);

    ##END CONNECTION
    disconnect(c);
    
  end

  disconnect(s);
endfunction


function msg = readstring (c)
  
  BUFFER_SIZE = 255;

  msg  = '';
  read = BUFFER_SIZE;

  while read >= BUFFER_SIZE
    newmsg = char(recv(c, BUFFER_SIZE));
    read = length(newmsg)
    msg = [ msg newmsg];
  end

endfunction

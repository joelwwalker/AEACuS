import logging

BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)

COLORS = {
    'WARNING'  : BLUE,
    'INFO'     : BLACK,
    'DEBUG'    : GREEN,
    'CRITICAL' : RED,
    'ERROR'    : RED,
    'RED'      : RED,
    'GREEN'    : GREEN,
    'YELLOW'   : YELLOW,
    'BLUE'     : BLUE,
    'MAGENTA'  : MAGENTA,
    'CYAN'     : CYAN,
    'WHITE'    : WHITE,
}

RESET_SEQ = "\033[0m"
COLOR_SEQ = "\033[1;%dm"
BOLD_SEQ  = "\033[1m"

class ColorFormatter(logging.Formatter):

    def __init__(self, *args, **kwargs):
        # can't do super(...) here because Formatter is an old school class
        logging.Formatter.__init__(self, *args, **kwargs)

    def format(self, record):
        levelname = record.levelname
        message   = logging.Formatter.format(self, record) + '$RESET'
        for k,v in COLORS.items():
            message = message.replace("$" + k,    COLOR_SEQ % (v+30))\
                         .replace("$BG" + k,  COLOR_SEQ % (v+40))\
                         .replace("$BG-" + k, COLOR_SEQ % (v+40))        
        
        
        if levelname == 'INFO':
            message   = message.replace("$RESET", '')\
                           .replace("$BOLD",  '')\
                           .replace("$COLOR", '')
            return message
        else:    
            color     = COLOR_SEQ % (30 + COLORS[levelname])
            message   = message.replace("$RESET", RESET_SEQ)\
                           .replace("$BOLD",  BOLD_SEQ)\
                           .replace("$COLOR", color)

        return message 

logging.ColorFormatter = ColorFormatter
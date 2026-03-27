from enum import IntEnum

class CountCategory(IntEnum):
    INTRON = 0
    UNSPLICED = 1
    SPLICED_UNANNOT = 2
    SPLICED_ANNOT = 3

class CountStrand(IntEnum):
    AMBIGUOUS = 0
    SENSE = 1
    ANTISENSE = 2

class CountType(IntEnum):
    """
    Defines designated integer categories for different types of counts.
    """
    INTRON_AMBIGUOUS = 0
    INTRON_SENSE = 1
    INTRON_ANTISENSE = 2
    UNSPLICED_AMBIGUOUS = 3
    UNSPLICED_SENSE = 4
    UNSPLICED_ANTISENSE = 5
    SPLICED_UNANNOT_AMBIGUOUS = 6
    SPLICED_UNANNOT_SENSE = 7
    SPLICED_UNANNOT_ANTISENSE = 8
    SPLICED_ANNOT_AMBIGUOUS = 9
    SPLICED_ANNOT_SENSE = 10
    SPLICED_ANNOT_ANTISENSE = 11
    
    # This private class attribute will cache the mapping.
    # It is initialized on the first call to from_name().
    __STR_TO_ENUM_MAP = None

    def to_column_name(self):
        """
        Converts the enum member to its corresponding column name.
        
        Example:
            CountType.INTRON_SENSE.to_column_name() -> 'intron_sense'
        """
        return self.name.lower()
    
    @classmethod
    def from_column_name(cls, name: str):
        """
        Gets an enum member from its lowercase string name.

        Example:
            CountType.from_name('spliced_sense') -> <CountType.SPLICED_SENSE: 7>
        """
        # The first time this method is called, it builds the mapping dictionary.
        if cls.__STR_TO_ENUM_MAP is None:
            cls.__STR_TO_ENUM_MAP = {e.to_column_name(): e for e in cls}       
        # All subsequent calls will use the cached dictionary for speed.
        return cls.__STR_TO_ENUM_MAP.get(name)

    @classmethod
    def columns(cls): 
        yield from (e.to_column_name() for e in cls)


! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

! $Revision$, $Date$
! $URL$
module ErrorMessages
  implicit none
  private
  ! Provides a uniform way to represent the success or failure of a calculation, 
  !   with optional text messages, as an alternative to, say, functions that return 
  !   integer codes with specific meanings.   
  
  ! Private (internal) parameters.  
  integer, parameter :: max_num_messages        = 100, &
                        max_message_length      = 256
                        
  ! -- Integer values that define the error state.
  !   The values are arbitrary and aren't available outside the module. 
  integer, parameter :: UndefinedState   = 0, &
                        SuccessState     = 1, &
                        WarningState     = 2, &
                        FailureState     = 3
  character (len = 21), parameter :: &
                        UndefinedMessageText = "Status is undefined."

  ! Public type: at least one variable of this type is declared in the 
  !   procedure that uses this module. The components of this type are 
  !   private, and so not accessable outside the module.  
  type, public :: ErrorMessage
    private
    
    ! lastMessage keeps track of the last message added to the pile
    ! currentMessage is used to iterate through a pile of messages
    integer                                :: lastMessage = 0, currentMessage = 0
    
    ! These members are allocated statically - this ensures that they'll 
    !   be available at run time and simplifies book keeping, but
    !   means that we can potentially run out of room.
    ! messageTexts keeps any textual information passed along with the status
    character (len = max_message_length),     &
             dimension(0:max_num_messages) :: messageTexts = trim(UndefinedMessageText)
             
    ! states holds one of the states we've defined
    integer, dimension(0:max_num_messages) :: states = UndefinedState
  end type ErrorMessage

  ! Here we declare which functions are visible outside the module. 
  !   This defines the interface. 
  ! -- Boolean functions for reporting the current error state
  public :: stateIsSuccess, &
            stateIsWarning, &
            stateIsFailure
 
  ! -- Routines for setting state and adding a new message
  public :: initializeState,   &
            setStateToSuccess, &
            setStateToWarning, &
            setStateToFailure, &
            setStateToCompleteSuccess
            
  ! -- Routines for iterating through the list of messages
  !  Pattern N ("Iterator") in the Gang of Four book. 
  public :: firstMessage,      &
            nextMessage,       &
            getCurrentMessage, & ! Gets the text associated with the state
            moreMessagesExist
  
  ! -- Introspection: inquire about the maximum lengths of strings
  public :: getErrorMessageLimits
  
  ! Usage: 
  
  ! program example
  !   type(ErrorMessage) :: status

  !   call doCalculation(..., status)
  !   if(.not. stateIsFailure(status)) call doAnotherCalculation(..., status)
  
  !   call firstMessage(status)
  !   historyLoop: do
  !     if (.not. moreMessagesExist(status)) exit historyLoop
  !     if(len_trim(getCurrentMessage(status)) > 0) print *, trim(getCurrentMessage(status))
  !     call nextMessage(status)
  !   end do historyLoop
  ! end program example

contains  
  ! -----------
  !  Functions for initializing or setting state with optional messages
  ! -----------
  function stateIsSuccess(messageVariable)
    type (ErrorMessage), intent ( in) :: messageVariable
    logical                           :: stateIsSuccess
    
    ! The current state is reported as success only if has been explicitly 
    !   set that way. No state may have ever been set (lastMessage = 0) 
    !   or the state might be something else (Warning, Failure). 
    stateIsSuccess =                        &
      messageVariable%lastMessage > 0 .and. &
      messageVariable%states(messageVariable%currentMessage) == SuccessState
  end function stateIsSuccess
  
  ! -----------
  function stateIsWarning(messageVariable)
    type (ErrorMessage), intent ( in) :: messageVariable
    logical                           :: stateIsWarning
    
    ! See notes in stateIsSuccess function
    stateIsWarning =                        &
      messageVariable%lastMessage > 0 .and. &
      messageVariable%states(messageVariable%currentMessage) == WarningState
  end function stateIsWarning
  
  ! -----------
  function stateIsFailure(messageVariable)
    type (ErrorMessage), intent ( in) :: messageVariable
    logical                           :: stateIsFailure
    
   ! See notes in stateIsSuccess function
    stateIsFailure =                        &
      messageVariable%lastMessage > 0 .and. &
      messageVariable%states(messageVariable%currentMessage) == FailureState
  end function stateIsFailure
  
  ! -----------
   function getCurrentMessage(messageVariable)
     type (ErrorMessage),    intent (in)  :: messageVariable
     character (len = max_message_length) :: getCurrentMessage

     ! If someone has never used one of the "set" functions they have
     !   no reason to be asking for the text. We give them some text 
     !   that tells them so. The state is undefined
     !   (not Success, Warning, or Failure). 
     if(messageVariable%currentMessage == 0) then
       getCurrentMessage = UndefinedMessageText
     else
       getCurrentMessage = trim(messageVariable%messageTexts(messageVariable%currentMessage))
     end if
   end function getCurrentMessage
   
  ! -----------
  !  Subroutines for initializing or setting state with optional messages
  ! -----------
  subroutine initializeState(messageVariable)
    type (ErrorMessage),           intent (out) :: messageVariable
    
    ! Starting with a clean slate. Everything is set to 
    !   undefined values
    messageVariable%currentMessage  = 0
    messageVariable%lastMessage     = 0
    messageVariable%messageTexts(:) = ""
    messageVariable%messageTexts(0) = trim(UndefinedMessageText)
    messageVariable%states(:)       = UndefinedState
  end subroutine initializeState
  
  ! -----------
  subroutine setStateToSuccess(messageVariable, messageText)
    type (ErrorMessage),           intent (inout) :: messageVariable
    character (len = *), optional, intent (   in) :: messageText
    
    ! Local variable
    integer :: thisMessage
    
    ! lastMessage keeps track of the end of the pile of messages; 
    !   currentMessage helps us iterate through the list. 
    ! We shouldn't change the state in the middle of iterating through
    !   the list. 
    thisMessage = messageVariable%lastMessage
    ! We only add a new message if we have enough room. We're most interested 
    !   in the first errors (presumably those lowest on the call stack) so we 
    !   overwrite the most recent messages if we don't have enough room. 
    if (thisMessage < max_num_messages) &
      thisMessage = thisMessage + 1
    messageVariable%currentMessage              = thisMessage
    ! When we set the state we go back to the top (most recent end)
    !   of the pile of messages. 
    messageVariable%lastMessage                 = thisMessage
    
    messageVariable%states(thisMessage)         = SuccessState
    if(present(messageText)) then
      messageVariable%messageTexts(thisMessage) = trim(messageText)
    else 
      messageVariable%messageTexts(thisMessage) = "" 
    end if
  end subroutine setStateToSuccess
  
  ! -----------
  subroutine setStateToWarning(messageVariable, messageText)
    type (ErrorMessage),           intent (inout) :: messageVariable
    character (len = *), optional, intent (   in) :: messageText
    
    ! Local variable
    integer :: thisMessage
    
    ! See comments in setStateToSuccess
    thisMessage = messageVariable%lastMessage
    if (thisMessage < max_num_messages) &
      thisMessage = thisMessage + 1
    messageVariable%lastMessage                 = thisMessage
    messageVariable%currentMessage              = thisMessage

    messageVariable%states(thisMessage)         = WarningState
    if(present(messageText)) then
      messageVariable%messageTexts(thisMessage) = trim(messageText)
    else 
      messageVariable%messageTexts(thisMessage) = "" 
    end if
  end subroutine setStateToWarning
  
  ! -----------
  subroutine setStateToFailure(messageVariable, messageText)
    type (ErrorMessage),           intent (inout) :: messageVariable
    character (len = *), optional, intent (   in) :: messageText
    
    ! Local variable
    integer :: thisMessage
    
    ! See comments in setStateToSuccess
    thisMessage = messageVariable%lastMessage
    if (thisMessage < max_num_messages) &
      thisMessage = thisMessage + 1
    messageVariable%lastMessage                 = thisMessage
    messageVariable%currentMessage              = thisMessage

    messageVariable%states(thisMessage)         = FailureState
    if(present(messageText)) then
      messageVariable%messageTexts(thisMessage) = trim(messageText)
    else 
      messageVariable%messageTexts(thisMessage) = "" 
    end if
  end subroutine setStateToFailure
  
  ! -----------
  subroutine setStateToCompleteSuccess(messageVariable, messageText)
    type (ErrorMessage),           intent (out) :: messageVariable
    character (len = *), optional, intent (in ) :: messageText

    ! Everything we've done is successful. We clear any messages 
    !   in the pile and tell everyone we're OK.     
    call initializeState(messageVariable)
    if(present(messageText)) then
      call setStateToSuccess(messageVariable, messageText) 
    else 
      call setStateToSuccess(messageVariable)
    end if
  end subroutine setStateToCompleteSuccess
  
  ! -----------
  !  Functions for iterating over a collection of messages
  ! -----------
   subroutine firstMessage(messageVariable)
     type (ErrorMessage), intent (inout) :: messageVariable
     
     messageVariable%currentMessage = min(1, messageVariable%lastMessage)
   end subroutine firstMessage
   
  ! -----------
   subroutine nextMessage(messageVariable)
     type (ErrorMessage), intent (inout) :: messageVariable
     
     ! Users should stop asking for the next message once we've reached the 
     !   top of the pile, but they might not. 
     ! if(messageVariable%currentMessage < messageVariable%lastMessage) &
       messageVariable%currentMessage = messageVariable%currentMessage + 1
   end subroutine nextMessage
   
  ! -----------
   function moreMessagesExist(messageVariable)
     type (ErrorMessage), intent (in) :: messageVariable
     logical                          :: moreMessagesExist
     
     ! Can we go on to the next message? 
     moreMessagesExist = &
       messageVariable%currentMessage <= messageVariable%lastMessage
   end function moreMessagesExist
   
  ! -----------
  !  Inquiry procedures
  ! -----------
  subroutine getErrorMessageLimits(messageVariable, maxNumberOfMessages, maxMessageLength)
     type (ErrorMessage), intent ( in) :: messageVariable
     integer, optional,   intent (out) :: maxNumberOfMessages, maxMessageLength
     
     ! Access to internal string lengths, maximum number of messages. 
     !  This is in the spirit of introspection, but isn't that useful since 
     !  everything is allocated at compile time. 
     if(present(maxNumberOfMessages)) &
       maxNumberOfMessages = max_num_messages
     if(present(maxMessageLength))    &
       maxMessageLength = max_message_length
  end subroutine getErrorMessageLimits
  
end module ErrorMessages

! -----------
! CVS/RCS modification history
!
! $Log: ErrorMessages.f95,v $
! Revision 1.7  2009/03/09 19:17:21  rpincus
! Placed code under GNU public license Version 2.
!
! Revision 1.6  2005/09/22 21:28:33  robert
! Corrected intent (changed to inout) in routines where old value is used.
!
! Revision 1.5  2005/09/10 03:57:45  robert
! Added CVS revision, date, and name tags so users can know which release
! they have.
!
! Revision 1.4  2005/02/12 00:43:26  robert
! Changed intents to inout in state setting functions - xlf is more
! strict than other compilers.
!
! Revision 1.3  2003/02/18 18:55:56  robert
! Changed name of finalization routine in monteCarloIllumination.
! Set status in read_Domain and getInfo_Domain in opticalProperties.
! Overwrote default message when setting state with no text in ErrorMessages.
!
! Revision 1.2  2003/02/17 19:17:56  robert
! Changed iteration scheme slightly; clarified comments.
!
! Revision 1.1.1.1  2003/02/14 00:18:10  robert
! I3RC Distribution
!
! -----------
